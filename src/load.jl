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

load_psp_file(dir::AbstractString, filename::AbstractString) = load_psp_file(joinpath(dir, filename))

"""
Load a pseudopotential file into its corresponding `AbstractPsP` subtype.
"""
load_psp(file::PsPFile) = formalism(file)(file)
load_psp(file::HghFile) = HghPsP(file)
load_psp(path::AbstractString) = return load_psp(load_psp_file(path))
load_psp(dir::AbstractString, filename::AbstractString) = load_psp(joinpath(dir, filename))

"""
List all known pseudopotential families.
"""
function list_families()
    artifacts_toml_path = Artifacts.find_artifacts_toml(pathof(PseudoPotentialIO))
    artifacts = Artifacts.load_artifacts_toml(artifacts_toml_path)
    artifact_names = keys(artifacts)
    return artifact_names
end

"""
List all known pseudopotential files in a given `family`, which can either be a path to 
a directory containing pseudopotential files or the name of a known pseudopotential family
artifact.
"""
function list_psp(family::AbstractString)
    if isdir(family)
        (_, _, files) = first(walkdir(family))
        psp_filenames = filter(f -> !startswith(f, "."), files)
        return psp_filenames
    else
        try
            family_path = @artifact_str(family)
            list_psp(family_path)
        catch
            error("PsP family $family does not exist")
        end
    end
end
