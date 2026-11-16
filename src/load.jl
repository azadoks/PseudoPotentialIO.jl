_FILE_EXT_LOADERS = Dict(".upf"  => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh"  => HghFile)

"""
Parse a pseudopotential file into a [`PsPFile`](@ref) struct.

Thread-safe: several tasks may read files at once. The error of a malformed UPF v2 file
read on a thread other than the first is only `unknown error of libxml2` (libxml2 prints
the details): EzXML installs its error handler on the thread that loads it.
"""
function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext in keys(_FILE_EXT_LOADERS) && return _FILE_EXT_LOADERS[ext](path)
    error("Unsupported PsP file extension $(ext)")
end
