function load_psp_file(path::AbstractString)
    _, ext = splitext(path)
    ext = lowercase(ext)
    ext == ".upf" && return UpfFile(path)
    ext == ".psp8" && return Psp8File(path)
    ext == ".hgh" && return HghFile(path)
    error("Unsupported file extension $(ext)")
end


function load_psp(path::AbstractString)
    file = load_psp_file(path)
    return load_psp(file)
end

function load_psp(file::PsPFile)
    formalism(file) == :norm_conserving && return NormConservingPsP(file)
    formalism(file) == :ultrasoft && return UltrasoftPsP(file)
    error("Unsupported formalism $(formalism(file))")
end

load_psp(file::HghFile) = HghPsP(file)
