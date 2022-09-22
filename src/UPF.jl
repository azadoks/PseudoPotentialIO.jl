module UPF

using EzXML

export UPF
include("UPF1.jl")

export UPF2
include("UPF2.jl")

export PSP8
include("PSP8.jl")

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
            return UPF1.parse_upf1(io)
        end
    elseif version == 2
        text = read(filename, String)
        # Remove end-of-file junk (input data, etc.)
        text = string(split(text, "</UPF>")[1], "</UPF>")
        # Clean any errant `&` characters
        text = replace(text, "&" => "")
        doc_root = root(parsexml(text))
        return UPF2.parse_upf2(doc_root)
    end
end

function get_psp_version(filename)
    version = 0
    open(filename, "r") do io
        for i = 1:2
            readline(io)
        end
        s = split(readline(io))
        version = parse(Int, s[1])
    end
    return version
end

function load_psp(filename::AbstractString)
    version = get_psp_version(filename)
    open(filename, "r") do io
        if version == 8
            psp = PSP8.parse_psp8(io)
            psp["header"]["filename"] = filename
            return psp
        else
            throw(ExceptionError("psp format version $(version) is not implemented."))
        end
    end
end

export load_upf
export load_psp

end
