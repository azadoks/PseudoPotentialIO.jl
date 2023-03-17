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

function load_psp_file(family_name_or_dir::AbstractString, filename::AbstractString)
    dir = _resolve_family(family_name_or_dir)
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
function load_psp(family_name_or_dir::AbstractString, filename::AbstractString)
    dir = _resolve_family(family_name_or_dir)
    path = joinpath(dir, filename)
    isfile(path) && return load_psp(path)
    return error("PsP $filename does not exist")
end

"""
List all known pseudopotential families with additional information about the number and
format(s) of the pseudopotentials within if requested (slow!).
"""
function list_families(;with_info=false)
    artifacts_toml_path = Artifacts.find_artifacts_toml(pathof(PseudoPotentialIO))
    artifacts = Artifacts.load_artifacts_toml(artifacts_toml_path)
    if with_info
        df = DataFrames.DataFrame(; name=[], n_psp=[], format=[])
        for (name, artifact) in artifacts
            if artifact_exists(Base.SHA1(artifact["git-tree-sha1"]))
                psps = load_family.(name)
                n_psp = size(psps, 1)
                format = unique(psps.Format)
                push!(df, (; name, n_psp, format))
            else
                push!(df, (; name, n_psp=missing, format=missing))
            end
        end
    else
        df = DataFrames.DataFrame(; name=collect(keys(artifacts)))
    end
    return df
end

"""
Show the elements contained in a pseudopotential family as a periodic table.
"""
function show_family_periodic_table(family_name_or_dir::AbstractString)
    family_df = load_family(family_name_or_dir)
    _show_family_header(family_df)
    show_family_periodic_table(family_df)
    return family_df
end
function show_family_periodic_table(family::DataFrames.DataFrame)
    elements = [PeriodicTable.elements[Symbol(element)] for element in family.Element]
    elements_struct = PeriodicTable.Elements(elements)

    println()
    show(Core.stdout, MIME("text/plain"), elements_struct)
    return nothing
end

"""
List the pseudopotentials in a pseudopotential family in a pretty table.
The elements for which pseudos are shown can be restricted by passing a list of strings,
e.g. ["Ag"].
"""
function show_family_list(family_name_or_dir::AbstractString; elements=[])
    family = load_family(family_name_or_dir)
    _show_family_header(family)
    show_family_list(family; elements)
    return family
end
function show_family_list(family::DataFrames.DataFrame; elements=[])
    if !isempty(elements)
        family = DataFrames.subset(family, :Element => e -> in.(e, Ref(elements)))
    end
    hlines = [0, 1]
    for i in 2:length(family.Element)
        if family.Element[i - 1] != family.Element[i]
            push!(hlines, i)
        end
    end
    push!(hlines, DataFrames.nrow(family) + 1)
    show(family; summary=false, allrows=true, eltypes=false, show_row_number=false,
         vlines=:all, hlines=hlines)
    println()
    return nothing
end

function _show_family_header(family::DataFrames.DataFrame)
    @printf "%10s: %s\n" "Family" family_name_or_dir
    @printf "%10s: %s\n" "Directory" _resolve_family(family_name_or_dir)
    @printf "%10s: %s\n" "Formalism" join(string.(unique(family.Formalism)), ", ")
    @printf "%10s: %d\n" "Pseudos" size(psp_df, 1)
    return nothing
end

"""
Load all pseudopotentials from a given family.
"""
function load_family(family_name_or_dir::AbstractString)
    dir = _resolve_family(family_name_or_dir)
    (_, _, filenames) = first(walkdir(dir))
    psp_filenames = filter(f -> !startswith(f, "."), filenames)
    psp_filenames = filter(f -> lowercase(splitext(f)[2]) in keys(_FILE_EXT_LOADERS),
                           psp_filenames)
    psp_paths = map(f -> joinpath(dir, f), psp_filenames)
    psp_files = load_psp_file.(psp_paths)
    df = DataFrames.DataFrame(; Element=element.(psp_files), Filename=psp_filenames,
                              Valence_Charge=valence_charge.(psp_files),
                              NLCC=has_core_density.(psp_files),
                              Spin_Orbit=has_spin_orbit.(psp_files),
                              Format=typeof.(psp_files),
                              Formalism=formalism.(psp_files),
                              File=psp_files)
    return df
end

# Resolve the name of a PsP family artifact into its directory
# or return the input if it's already a directory
function _resolve_family(family_name_or_dir::AbstractString)
    isdir(family_name_or_dir) && return family_name_or_dir
    try
        dir = @artifact_str(family_name_or_dir)
        return dir
    catch
        error("PsP family $family does not exist")
    end
end
