module UPF

using EzXML

include("upf1.jl")
export parse_upf1

include("upf2.jl")
export parse_upf2

function get_upf_version(filename)
    open(filename, "r") do io
        line = readline(io)
        if occursin("<PP_INFO>", line)
            return 1
        elseif occursin("UPF version", line)
            return 2
        else
            throw(Error("Unknown UPF version"))
        end
    end
end

function load_upf(filename::AbstractString)
    version = get_upf_version(filename)
    if version == 1
        open(filename, "r") do io
            return parse_upf1(io)
        end
    elseif version == 2
        open(filename, "r") do io
            text = read(io, String)
        end
        # Remove end-of-file junk
        text = string(split(text, "</UPF>")[1], "</UPF>")
        # Clean any errant "&" characters
        text = replace(text, "&" => "")
        doc_root = root(parsexml(text))
        return parse_upf2(doc_root)
    end
end

export load_upf

end
