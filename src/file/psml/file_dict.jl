#= PSML <-> nested dictionaries

The same layout as `upf/file_dict_v2.jl`, keyed by the tag names as they are (`"pseudo-atom-spec"`):
- an element with only numbers (`data`, `grid-data`) is a `Vector{Float64}`
- an element with attributes or children is an `OrderedDict`: attributes are scalar entries
  (`Int`, `Float64`, `Bool` or `String`), children are array or dictionary entries, and the
  element's text (if any, e.g. `input-file`) is under "content"
- elements that the schema lets repeat (`provenance`, `shell`, `slps`, `proj`, ...) are
  always a `Vector` of dictionaries, in file order

The sizes (`npts` of `grid` and `data`) are checked against the data, dropped, and written
back from the data. The order of the entries is the order of the file, and the writer writes
them in the order of the dictionary.

Based on the schema (libpsml `doc/schema/psml.rnc`) and libpsml's parser
(`m_psml_parsing_helpers.F90`):
  - `yes`/`no` attributes are `Bool`s; `annotation` attributes are free text
  - a `data` uses the grid of its `radfunc`, else of the closest enclosing element, else the
    top-level one; it has the grid's size unless its `npts` says it is shorter
Nothing else is enforced (element order, required elements, version), so older files
(PSML 0.8 `header`, `pseudopotential-operator`, `projectors`) read as well.
=#

const PsmlDict = OrderedDict{String,Any}

# elements that can repeat in their parent
const PSML_LIST_ELEMENTS = Set(["provenance", "input-file", "shell", "functional",
                                "semilocal-potentials", "slps", "nonlocal-projectors",
                                "projectors", "proj", "pseudo-wave-functions", "pswf"])
# attributes that are text even when they look like numbers
const PSML_TEXT_ATTRIBUTES = Set(["version", "uuid", "xmlns", "energy_unit", "length_unit",
                                  "creator", "date", "name", "atomic-label", "relativity",
                                  "flavor", "type", "set", "l"])
# attributes that are booleans
const PSML_BOOL_ATTRIBUTES = Set(["core-corrections", "spin-dft", "polarized", "meta-gga",
                                  "is-unscreening-charge", "rescaled-to-z-pseudo",
                                  "is-unscreening-tau"])

function parse_psml_attribute(tag, key, value::AbstractString)
    (tag == "annotation" || key in PSML_TEXT_ATTRIBUTES) && return String(value)
    v = strip(value)
    key in PSML_BOOL_ATTRIBUTES && v in ("yes", "no") && return v == "yes"
    x = tryparse(Int, v)
    !isnothing(x) && return x
    x = parse_fortran_reals(v)
    return isnothing(x) || length(x) != 1 ? String(value) : only(x)
end

format_psml_attribute(x::Bool) = x ? "yes" : "no"
format_psml_attribute(x) = string(x)

"""
Read a PSML file (its text) into nested dictionaries (see `src/file/psml/file_dict.jl`).
"""
function psml_to_dict(text::AbstractString)::PsmlDict
    root = only(n for n in XML.children(XML.parse(XML.Node, text))
                if XML.nodetype(n) == XML.Element)
    XML.tag(root) == "psml" || error("Not a PSML file: the root element is $(XML.tag(root))")
    return psml_element_to_dict(root, nothing)
end

psml_elements(node) = [c for c in XML.children(node) if XML.nodetype(c) == XML.Element]
psml_attributes(node) = something(XML.attributes(node), Dict{String,String}())

function grid_npts(node)
    npts = tryparse(Int, get(psml_attributes(node), "npts", ""))
    (isnothing(npts) || npts <= 0) && error("PSML <grid> without a valid npts")
    return npts
end

# `grid`: the number of points of the enclosing grid (`nothing` if there is none yet)
function psml_element_to_dict(node, grid)
    tag = XML.tag(node)
    attributes = psml_attributes(node)
    children = psml_elements(node)
    text = join(XML.value(c) for c in XML.children(node)
                if XML.nodetype(c) in (XML.Text, XML.CData))
    if tag in ("data", "grid-data")
        values = parse_xml_numbers(text)
        isnothing(values) && error("PSML <$tag> has non-numeric values")
        isempty(values) && error("PSML <$tag> is empty")
        isnothing(grid) && error("PSML <$tag> without a grid")
        npts = haskey(attributes, "npts") ? parse(Int, attributes["npts"]) : grid
        tag == "data" && npts > grid &&
            error("PSML <data> has npts=$npts, but its grid has $grid points")
        length(values) == npts ||
            error("PSML <$tag> has $(length(values)) values, but should have $npts")
        return values
    end

    d = PsmlDict()
    for (k, v) in attributes
        tag == "grid" && k == "npts" && continue  # it is the length of grid-data
        d[k] = parse_psml_attribute(tag, k, v)
    end
    isempty(strip(text)) || (d["content"] = text)
    if tag == "grid"
        grid = grid_npts(node)
    else  # a grid of this element applies to its other children
        for c in children
            XML.tag(c) == "grid" && (grid = grid_npts(c))
        end
    end
    for child in children
        key = XML.tag(child)
        value = psml_element_to_dict(child, grid)
        if key in PSML_LIST_ELEMENTS
            push!(get!(d, key, Any[]), value)
        elseif haskey(d, key)  # repeated element the schema does not list
            d[key] isa Vector{Any} || (d[key] = Any[d[key]])
            push!(d[key], value)
        else
            d[key] = value
        end
    end
    return d
end

function psml_dict_to_element(tag::AbstractString, value)
    if is_numeric_array(value)
        attributes = tag == "data" ? (; npts=string(length(value))) : (;)
        return XML.Element(tag, XML.Text(format_xml_numbers(value)); attributes...)
    end
    attributes = Pair{String,String}[]
    if tag == "grid" && haskey(value, "grid-data")
        push!(attributes, "npts" => string(length(value["grid-data"])))
    end
    children = []
    for (k, v) in value
        if isnothing(v)
            continue
        elseif k == "content"
            node = occursin("]]>", v) ? XML.Text(v) : XML.CData(v)
            push!(children, node)
        elseif v isa AbstractDict || is_numeric_array(v)
            push!(children, psml_dict_to_element(k, v))
        elseif v isa AbstractVector
            append!(children, psml_dict_to_element.(k, v))
        else
            push!(attributes, k => format_psml_attribute(v))
        end
    end
    return XML.Element(tag, children...; (Symbol(k) => v for (k, v) in attributes)...)
end

"""
Write nested dictionaries (see `psml_to_dict`) as a PSML file to `io`.
"""
function dict_to_psml(io::IO, d::AbstractDict)
    println(io, """<?xml version="1.0" encoding="UTF-8"?>""")
    XML.write(io, psml_dict_to_element("psml", d))
    println(io)
end
