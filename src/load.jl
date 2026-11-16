_FILE_EXT_LOADERS = Dict(".upf"  => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh"  => HghFile,
                         ".psml" => PsmlFile)

"""
Parse a pseudopotential file into a [`PsPFile`](@ref) struct. A value that does not follow
the specification of the format is an error, a warning, or nothing, by `validation`
(`:error`, `:warn` or `:none`; see [`VALIDATION`](@ref)).

Thread-safe: several tasks may read files at once.
"""
function load_psp_file(path::AbstractString; validation=VALIDATION[])
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext in keys(_FILE_EXT_LOADERS) && return _FILE_EXT_LOADERS[ext](path; validation)
    error("Unsupported PsP file extension $(ext)")
end
