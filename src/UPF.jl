module UPF

using EzXML

export UPF
include("UPF1.jl")

export UPF2
include("UPF2.jl")

export PSP8
include("PSP8.jl")

"""
    get_upf_version(filename::AbstractString)

Parse the version of a UPF file from its contents.
Supported versions are old UPF (v1) and UPF with schema (v2).

!!! Note
UPF v2 without schema is not supported.
"""
function get_upf_version(filename::AbstractString)::Int
    open(filename, "r") do io
        line = readline(io)
        # Old UPF files start with the `<PP_INFO>` section
        if occursin("<PP_INFO>", line)
            return 1
        # New UPF files with schema are in XML and start with a version tag
        elseif occursin("UPF version=\"2.0.1\"", line)
            return 2
        else
            throw(Error("Unknown UPF version"))
        end
    end
end

"""
    load_upf(filename::AbstractString)

Load a UPF file, supports old (v1) and v2 with schema.
"""
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

"""
    get_psp_version(filename::AbstractString)

Parse the version of a PSP file from its contents.
"""
function get_psp_version(filename::AbstractString)::Int
    open(filename, "r") do io
        # Discard the first two lines
        for _ = 1:2
            readline(io)
        end
        # The version number is the first entry on line 3
        s = split(readline(io))
        return parse(Int, s[1])
    end
end

"""
    load_psp(filename::AbstractString)

Load a PSP file, supports only version 8.
"""
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
