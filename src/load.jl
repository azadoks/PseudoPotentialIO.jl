#TODO functions for listing pseudo families, pseudos of a given element, etc.
#TODO functions for working with the pseudo family artifacts

"""
Parse a pseudopotential file into a `PsPFile` struct. 
"""
function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext == ".upf" && return UpfFile(path)
    ext == ".psp8" && return Psp8File(path)
    ext == ".hgh" && return HghFile(path)
    error("Unsupported file extension $(ext)")
end

"""
Parse pseudopotential file contents into an `AbstractPsP` struct.
"""
function load_psp(file::PsPFile)
    formalism(file) == :norm_conserving && return NormConservingPsP(file)
    formalism(file) == :ultrasoft && return UltrasoftPsP(file)
    error("Unsupported formalism $(formalism(file))")
end
load_psp(file::HghFile) = HghPsP(file)

load_psp(path::AbstractString) = return load_psp(load_psp_file(path))
load_psp(dir::AbstractString, filename::AbstractString) = load_psp(joinpath(dir, filename))
