#= UPF v2 <-> nested dictionaries

Every element becomes an entry of its parent, keyed by its tag in lower case without `PP_`
and without the index suffix (`PP_BETA.2` -> "beta"):
- an element with only numbers is a `Vector{Float64}`
- an element with attributes or children is an `OrderedDict`: attributes are scalar entries
  (`Int`, `Float64`, `Bool` or `String`), children are array or dictionary entries, and the
  element's text (if any) is under "content"
- indexed elements (`PP_CHI.1`, `PP_QIJ.1.2`, ...) are a `Vector` of dictionaries, in file
  order; the suffix is stored in their attributes (`index`, or `first_index`, ...), which
  is how the writer names them again

The bookkeeping attributes of arrays (`size`, `type`, `columns`, `rows`) are checked against
the data, dropped, and written back from the data. The order of the entries is the order
of the file, and the writer writes them in the order of the dictionary.

Additional behavior based on QE's upflib:
  - elements are looked up by the name with the suffix, which wins over an `index` attribute
    ("absent or incorrect in some UPF v.2 files")
  - logicals are read list-directed (`T`, `.TRUE.`, `true`,...)
  - numbers are list-directed (`1.0D-3`, `2.0-319`, `3*0.0`).
=#

const UpfDict = OrderedDict{String,Any}

# attributes that are text even when they look like numbers or logicals (`element="F"`,
# `date="180509"`)
const UPF_TEXT_ATTRIBUTES = Set(["generated", "author", "date", "comment", "element",
                                 "pseudo_type", "relativistic", "functional", "label",
                                 "els", "shape", "breakpoints", "breakpoint_indices"])
# attributes that only describe the layout of an array
const UPF_ARRAY_ATTRIBUTES = ("size", "type", "columns", "rows")
# the attributes that the index suffix of an element is (default: "index")
const UPF_SUFFIX_ATTRIBUTES = Dict("qij" => ("first_index", "second_index"),
                                   "qijl" => ("first_index", "second_index",
                                              "angular_momentum"))
suffix_attributes(key) = get(UPF_SUFFIX_ATTRIBUTES, key, ("index",))

function parse_upf_attribute(key::AbstractString, value::AbstractString)
    value = strip(value)
    key in UPF_TEXT_ATTRIBUTES && return String(value)
    occursin(r"^\.?(T|F|TRUE|FALSE)\.?$"i, value) && return uppercase(lstrip(value, '.'))[1] == 'T'
    x = tryparse(Int, value)
    !isnothing(x) && return x
    x = parse_fortran_reals(value)
    return isnothing(x) || length(x) != 1 ? String(value) : only(x)
end

"""
Read a UPF v2 file (its text) into nested dictionaries (see `src/file/upf/file_dict_v2.jl`).
"""
function upf_to_dict(text::AbstractString)::UpfDict
    # junk after the root (input data, ...) and unescaped ampersands are common
    text = escape_ampersand(first(split(text, "</UPF>")) * "</UPF>")
    root = only(n for n in XML.children(XML.parse(XML.Node, text))
                if XML.nodetype(n) == XML.Element)
    XML.tag(root) == "UPF" || error("Not a UPF v2 file: the root element is $(XML.tag(root))")
    return upf_element_to_dict(root)
end

function upf_element_to_dict(node)
    d = UpfDict()
    attributes = something(XML.attributes(node), Dict{String,String}())
    for (k, v) in attributes
        k in UPF_ARRAY_ATTRIBUTES || (d[k] = parse_upf_attribute(k, v))
    end
    children = [c for c in XML.children(node) if XML.nodetype(c) == XML.Element]
    text = join(XML.value(c) for c in XML.children(node)
                if XML.nodetype(c) in (XML.Text, XML.CData))
    if !isempty(strip(text)) || isempty(children)
        numbers = parse_xml_numbers(text)
        content = isnothing(numbers) ? text : numbers
        if haskey(attributes, "size") && content isa Vector &&
           length(content) != parse(Int, attributes["size"])
            error("UPF $(XML.tag(node)) has $(length(content)) values, but size=$(attributes["size"])")
        end
        # text keeps its dictionary so that scalars are always attributes
        content isa Vector && isempty(d) && isempty(children) && return content
        isempty(content) || (d["content"] = content)
    end
    for child in children
        base, suffix... = split(XML.tag(child), '.')
        key = lowercase(replace(base, r"^PP_" => ""))
        value = upf_element_to_dict(child)
        indices = tryparse.(Int, suffix)
        if !isempty(indices) && all(!isnothing, indices)
            value isa AbstractDict || (value = UpfDict("content" => value))
            for (attribute, index) in zip(suffix_attributes(key), indices)
                value[attribute] = index  # the suffix is what QE goes by
            end
            push!(get!(d, key, Any[]), value)
        elseif haskey(d, key)  # repeated element without an index
            d[key] isa Vector{Any} || (d[key] = Any[d[key]])
            push!(d[key], value)
        else
            d[key] = value
        end
    end
    return d
end

format_upf_attribute(x::Bool) = x ? "T" : "F"
format_upf_attribute(x) = string(x)  # shortest text that reads back to the same value

function upf_dict_to_element(tag::AbstractString, value)
    if is_numeric_array(value)
        return XML.Element(
            tag, XML.Text(format_xml_numbers(value)); size=string(length(value))
        )
    end
    attributes = Pair{String,String}[]
    children = []
    for (k, v) in value
        if isnothing(v)
            continue
        elseif k == "content"
            is_numeric_array(v) && push!(attributes, "size" => string(length(v)))
            push!(children, XML.Text(is_numeric_array(v) ? format_xml_numbers(v) : v))
        elseif v isa AbstractDict || is_numeric_array(v)
            push!(children, upf_dict_to_element("PP_" * uppercase(k), v))
        elseif v isa AbstractVector
            for item in v
                name = "PP_" * uppercase(k)
                if item isa AbstractDict
                    indices = [get(item, a, nothing) for a in suffix_attributes(k)]
                    all(!isnothing, indices) && (name *= "." * join(indices, '.'))
                end
                push!(children, upf_dict_to_element(name, item))
            end
        else
            push!(attributes, k => format_upf_attribute(v))
        end
    end
    return XML.Element(tag, children...; (Symbol(k) => v for (k, v) in attributes)...)
end

"""
Write nested dictionaries (see `upf_to_dict`) as a UPF v2 file to `io`.
"""
function dict_to_upf(io::IO, d::AbstractDict)
    println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
    XML.write(io, upf_dict_to_element("UPF", d))
    println(io)
end
