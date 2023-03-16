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
        df = DataFrames.DataFrame(; name=[], n_psp=[], format=[])
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

function show_family_table(family_or_dir::AbstractString)
    psp_df = _show_family_header(family_or_dir)
    elements = [PeriodicTable.elements[Symbol(element)] for element in psp_df.Element]
    elements_struct = PeriodicTable.Elements(elements)

    println()
    show(Core.stdout, MIME("text/plain"), elements_struct)
    return nothing
end

function show_family_list(family_or_dir::AbstractString; elements=[])
    psp_df = _show_family_header(family_or_dir)
    if !isempty(elements)
        psp_df = DataFrames.subset(psp_df, :Element => e -> in.(e, Ref(elements)))
    end
    hlines = [0, 1]
    for i in 2:length(psp_df.Element)
        if psp_df.Element[i-1] != psp_df.Element[i]
            push!(hlines, i)
        end
    end
    push!(hlines, DataFrames.nrow(psp_df) + 1)
    show(psp_df; summary=false, allrows=true, eltypes=false, show_row_number=false,
         vlines=:all, hlines=hlines)
    println()
    return nothing
end

function _show_family_header(family_or_dir::AbstractString)
    psp_df = list_psp(family_or_dir)
    @printf "%10s: %s\n" "Family" family_or_dir
    @printf "%10s: %s\n" "Directory" _resolve_family(family_or_dir)
    @printf "%10s: %s\n" "Formalism" join(string.(unique(psp_df.Formalism)), ", ") 
    @printf "%10s: %d\n" "Pseudos" size(psp_df, 1)
    return psp_df
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
    df = DataFrame(; Element=element.(psp_files), Filename=psp_filenames,
                   Valence_Charge=valence_charge.(psp_files),
                   NLCC=has_core_density.(psp_files),
                   Spin_Orbit=has_spin_orbit.(psp_files),
                   Format=typeof.(psp_files),
                   Formalism=formalism.(psp_files))
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
