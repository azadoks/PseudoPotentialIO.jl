_FILE_EXT_LOADERS = Dict(".upf" => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh" => HghFile)

"""
Parse a pseudopotential file into a `PsPFile` struct. 
"""
function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext in keys(_FILE_EXT_LOADERS) && return _FILE_EXT_LOADERS[ext](path)
    return error("Unsupported PsP file extension $(ext)")
end

function load_psp_file(family_or_dir::AbstractString, filename::AbstractString)
    dir = _resolve_family(family_or_dir)
    path = joinpath(dir, filename)
    isfile(path) && return load_psp_file(path)
    return error("PsP $filename does not exist")
end

"""
Load a pseudopotential file into its corresponding `AbstractPsP` subtype.
"""
load_psp(file::PsPFile) = formalism(file)(file)
load_psp(file::HghFile) = HghPsP(file)
load_psp(path::AbstractString) = return load_psp(load_psp_file(path))
function load_psp(family_or_dir::AbstractString, filename::AbstractString)
    dir = _resolve_family(family_or_dir)
    path = joinpath(dir, filename)
    isfile(path) && return load_psp(path)
    return error("PsP $filename does not exist")
end

"""
List all known pseudopotential families with additional information about the number and
format(s) of the pseudopotentials within if requested (slow!).
"""
function list_families(with_info=false)
    artifacts_toml_path = Artifacts.find_artifacts_toml(pathof(PseudoPotentialIO))
    artifacts = Artifacts.load_artifacts_toml(artifacts_toml_path)
    if with_info
        df = DataFrame(; name=[], n_psp=[], format=[])
        for (name, artifact) in artifacts
            if artifact_exists(Base.SHA1(artifact["git-tree-sha1"]))
                psps = list_psp.(name)
                n_psp = size(psps, 1)
                format = unique(psps.format)
                push!(df, (; name, n_psp, format))
            else
                push!(df, (; name, n_psp=missing, format=missing))
            end
        end
    else
        df = DataFrame(; name=collect(keys(artifacts)))
    end
    return df
end

"""
List all known pseudopotential files in a given `family`, which can either be a path to 
a directory containing pseudopotential files or the name of a known pseudopotential family
artifact.
"""
function list_psp(family_or_dir::AbstractString)
    dir = _resolve_family(family_or_dir)
    (_, _, filenames) = first(walkdir(dir))
    psp_filenames = filter(f -> !startswith(f, "."), filenames)
    psp_filenames = filter(f -> lowercase(splitext(f)[2]) in keys(_FILE_EXT_LOADERS),
                           psp_filenames)
    psp_paths = map(f -> joinpath(dir, f), psp_filenames)
    psp_files = load_psp_file.(psp_paths)
    df = DataFrame(; filename=psp_filenames, element=element.(psp_files),
                   valence_charge=valence_charge.(psp_files),
                   NLCC=has_core_density.(psp_files),
                   spin_orbit=has_spin_orbit.(psp_files), format=typeof.(psp_files))
    return df
end

# Resolve the name of a PsP family artifact into its directory
# or return the input if it's already a directory
function _resolve_family(family_or_dir::AbstractString)
    isdir(family_or_dir) && return family_or_dir
    try
        dir = @artifact_str(family_or_dir)
        return dir
    catch
        error("PsP family $family does not exist")
    end
end