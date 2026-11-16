#= UPF <-> `UpfFile`, through the dictionaries of `file_dict_v2.jl` (UPF v1 is read into
them too): `_from_dict` builds the structs (and checks what the format leaves to the
reader), `_to_dict` is its inverse. Most structs are their element's attributes, field for
field (`upf_struct`, `upf_attributes`). =#

upf2_parse_psp(io::IO; identifier="") =
    _from_dict(UpfFile, upf_to_dict(read(io, String)); identifier)

upf2_dump_psp(io::IO, file::UpfFile) = dict_to_upf(io, _to_dict(file))

# Entry `key` of `d` as a `T`. A missing logical is false, as in QE; anything else missing
# is `nothing` if `T` allows it and an error if not.
function upf_get(::Type{T}, d, key::AbstractString) where {T}
    S = Base.nonnothingtype(T)
    x = get(d, key, nothing)
    x isa AbstractDict && S <: AbstractVector && (x = get(x, "content", nothing))
    if isnothing(x)
        T === Bool && return false
        Nothing <: T && return nothing
        error("UPF: `$key` is missing")
    end
    try
        return S === String ? string(x) : convert(S, x)  # e.g. `l="2.0"` is an integer
    catch
        error("UPF: cannot read `$key` = $(repr(x)) as $S")
    end
end

# A struct from its element's attributes, by field name; `fields` gives other values.
upf_struct(::Type{T}, d; fields...) where {T} =
    T((haskey(fields, f) ? fields[f] : upf_get(fieldtype(T, f), d, string(f))
       for f in fieldnames(T))...)

# The attributes of an element from a struct, and its array as the content.
function upf_attributes(x, array=nothing)
    d = UpfDict(string(f) => getfield(x, f) for f in fieldnames(typeof(x))
                if f != array && !isnothing(getfield(x, f)))
    isnothing(array) || (d["content"] = getfield(x, array))
    return d
end

# the elements `key` of `d` (indexed or not) as a vector; `d` may be an empty element
function upf_list(d, key)
    d isa AbstractDict && haskey(d, key) || return []
    x = d[key]
    return x isa AbstractVector && !is_numeric_array(x) ? x : [x]
end

# an element as a dictionary, also when it is empty or missing
upf_element(d, key) = (x = get(d, key, nothing); x isa AbstractDict ? x : UpfDict())

upf_array(d, key) = upf_get(Union{Nothing,Vector{Float64}}, d, key)
nonempty(x) = isnothing(x) || isempty(x) ? nothing : x

# The names of the tables that can have breakpoints, with their dictionaries (see
# `UpfFile`): the arrays of `PP_UPF`, of `PP_NONLOCAL` and of `PP_PSWFC`.
function upf_tables(d)
    tables = Pair{String,Any}[]
    for parent in (d, get(d, "nonlocal", nothing), get(d, "pswfc", nothing))
        parent isa AbstractDict || continue
        for (key, value) in parent
            for (i, item) in enumerate(upf_list(parent, key))
                (is_numeric_array(item) || item isa AbstractDict && haskey(item, "content")) ||
                    continue
                suffix = item isa AbstractDict && haskey(item, "index") ? ".$(item["index"])" : ""
                push!(tables, "PP_$(uppercase(key))$suffix" => (parent, key, i))
            end
        end
    end
    return tables
end

parse_radii(text) = parse.(Float64, split(replace(uppercase(text), "D" => "E")))

function _from_dict(::Type{UpfFile}, d::AbstractDict; identifier="")
    header = _from_dict(UpfHeader, d["header"])
    header.pseudo_type == "SL" && error("Semilocal pseudopotentials are not supported")
    mesh = _from_dict(UpfMesh, d["mesh"])
    local_ = header.is_coulomb ? nothing : upf_array(d, "local")
    rhoatom = upf_get(Vector{Float64}, d, "rhoatom")
    nonlocal = _from_dict(UpfNonlocal, upf_element(d, "nonlocal"))
    # projectors that are zero are dropped (see `_from_dict(UpfNonlocal)`): so are they from
    # the count (the rest of the header is validated already)
    nproj = length(nonlocal.betas)
    if nproj < header.number_of_proj
        fields = (f === :number_of_proj ? nproj : getfield(header, f)
                  for f in fieldnames(UpfHeader))
        header = with_validation(() -> UpfHeader(fields...), :none)
    end
    full_wfc = haskey(d, "full_wfc") ? _from_dict(UpfFullWfc, d["full_wfc"]) : nothing
    spin_orb = haskey(d, "spin_orb") ? _from_dict(UpfSpinOrb, d["spin_orb"]) : nothing
    paw = haskey(d, "paw") ? _from_dict(UpfPaw, d["paw"]) : nothing
    gipaw = haskey(d, "gipaw") ? _from_dict(UpfGipaw, d["gipaw"]) : nothing
    breakpoints = Dict{String,Vector{Float64}}()
    for (name, (parent, key, i)) in upf_tables(d)
        item = upf_list(parent, key)[i]
        item isa AbstractDict && haskey(item, "breakpoints") &&
            (breakpoints[name] = parse_radii(item["breakpoints"]))
    end
    return UpfFile(identifier, string(d["version"]), upf_info(get(d, "info", nothing)), header,
                   mesh, nonempty(upf_array(d, "nlcc")), upf_array(d, "taumod"), local_,
                   nonlocal,
                   [_from_dict(UpfChi, c) for c in upf_list(get(d, "pswfc", nothing), "chi")],
                   full_wfc, rhoatom, upf_array(d, "tauatom"), spin_orb, paw, gipaw,
                   breakpoints)
end

# the text of PP_INFO, with that of PP_INPUTFILE
upf_info(::Nothing) = nothing
upf_info(::AbstractVector) = ""  # empty
upf_info(info::AbstractDict) =
    get(info, "content", "") * upf_info(get(info, "inputfile", Float64[]))

function _to_dict(file::UpfFile)
    d = UpfDict("version" => "2.0.1")
    isnothing(file.info) || (d["info"] = UpfDict("content" => file.info))
    d["header"] = upf_attributes(file.header)
    d["mesh"] = _to_dict(file.mesh)
    d["nlcc"] = file.nlcc
    d["local"] = file.local_
    d["nonlocal"] = _to_dict(file.nonlocal)
    d["pswfc"] = UpfDict("chi" => [upf_attributes(c, :chi) for c in file.pswfc])
    isnothing(file.full_wfc) || (d["full_wfc"] = _to_dict(file.full_wfc))
    d["rhoatom"] = file.rhoatom
    d["taumod"] = file.taumod
    d["tauatom"] = file.tauatom
    isnothing(file.spin_orb) || (d["spin_orb"] = _to_dict(file.spin_orb))
    isnothing(file.paw) || (d["paw"] = _to_dict(file.paw))
    isnothing(file.gipaw) || (d["gipaw"] = _to_dict(file.gipaw))
    # breakpoints (an extension, see `UpfFile`), with the number of mesh points below each
    tables = Dict(upf_tables(d))
    for (name, breaks) in file.breakpoints
        haskey(tables, name) ||
            throw(ArgumentError("breakpoints of $name, which the file does not have"))
        parent, key, i = tables[name]
        item = upf_list(parent, key)[i]
        item isa AbstractDict || (item = parent[key] = UpfDict("content" => item))
        item["breakpoints"] = join((@sprintf("%.16E", b) for b in breaks), ' ')
        item["breakpoint_indices"] = join((searchsortedfirst(file.mesh.r, b) - 1
                                           for b in breaks), ' ')
    end
    return d
end

function _from_dict(::Type{UpfHeader}, d)
    # QE writes the short and the long names separated by spaces or dashes ("SLA-PW-PBX-PBC")
    functional = strip(d["functional"])
    startswith(uppercase(functional), "XC-") ||
        (functional = join(split(functional, x -> isspace(x) || x == '-'; keepempty=false), ' '))
    return upf_struct(UpfHeader, d; functional)
end

function _from_dict(::Type{UpfMesh}, d)
    r = upf_get(Vector{Float64}, d, "r")
    rab = upf_get(Vector{Float64}, d, "rab")
    mesh = something(upf_get(Union{Nothing,Int}, d, "mesh"), length(r))
    return upf_struct(UpfMesh, d; r, rab, mesh)
end

function _to_dict(mesh::UpfMesh)
    d = upf_attributes(mesh)
    d["r"] = pop!(d, "r")
    d["rab"] = pop!(d, "rab")
    return d
end

function _from_dict(::Type{UpfNonlocal}, d)
    betas = [upf_struct(UpfBeta, b; beta=upf_get(Vector{Float64}, b, "content"))
             for b in upf_list(d, "beta")]
    nbeta = length(betas)
    dij = isempty(betas) ? zeros(0, 0) : upf_get(Vector{Float64}, d, "dij")
    length(dij) == nbeta^2 ||
        error("UPF: PP_DIJ has $(length(dij)) values for $nbeta projectors")
    dij = reshape(dij, nbeta, nbeta)
    # projectors that are zero are dropped
    keep = [!all(iszero, b.beta) for b in betas]
    augmentation = haskey(d, "augmentation") ?
                   _from_dict(UpfAugmentation, d["augmentation"], nbeta) : nothing
    return UpfNonlocal(betas[keep], dij[keep, keep], augmentation)
end

function _to_dict(nl::UpfNonlocal)
    d = UpfDict("beta" => [upf_attributes(b, :beta) for b in nl.betas], "dij" => vec(nl.dij))
    isnothing(nl.augmentation) || (d["augmentation"] = _to_dict(nl.augmentation))
    return d
end

function _from_dict(::Type{UpfAugmentation}, d, nbeta)
    q = upf_get(Vector{Float64}, d, "q")
    length(q) == nbeta^2 || error("UPF: PP_Q has $(length(q)) values for $nbeta projectors")
    q = reshape(q, nbeta, nbeta)
    # PP_QFCOEF is the flat Fortran array qfcoef(nqf, nqlc, nbeta, nbeta), symmetric in the
    # projector indices; stored as one `UpfQfcoef` (nqf × nqlc, flattened) per pair i ≤ j,
    # as in UPF v1
    qfcoefs = upf_array(d, "qfcoef")
    if !isnothing(qfcoefs)
        qfcoef = reshape(qfcoefs, upf_get(Int, d, "nqf"), :, nbeta, nbeta)
        qfcoefs = [UpfQfcoef(vec(qfcoef[:, :, i, j]), nothing, i, j)
                   for i in 1:nbeta for j in i:nbeta]
    end
    qijs = [upf_struct(UpfQij, x; qij=upf_get(Vector{Float64}, x, "content"))
            for x in upf_list(d, "qij")]
    qijls = [upf_struct(UpfQijl, x; qijl=upf_get(Vector{Float64}, x, "content"))
             for x in upf_list(d, "qijl")]
    return upf_struct(UpfAugmentation, d; q, qfcoefs, multipoles=upf_array(d, "multipoles"),
                      rinner=upf_array(d, "rinner"), qijs=nonempty(qijs),
                      qijls=nonempty(qijls))
end

function _to_dict(aug::UpfAugmentation)
    d = upf_attributes(aug)
    for k in ("q", "multipoles", "qfcoefs", "rinner", "qijs", "qijls")
        delete!(d, k)
    end
    d["q"] = vec(aug.q)
    d["multipoles"] = aug.multipoles
    if !isnothing(aug.qfcoefs)  # the flat Fortran array, see `_from_dict`
        nq = size(aug.q, 1)
        qfcoef = zeros(length(first(aug.qfcoefs).qfcoef), nq, nq)
        for c in aug.qfcoefs
            qfcoef[:, c.first_index, c.second_index] = c.qfcoef
            qfcoef[:, c.second_index, c.first_index] = c.qfcoef
        end
        d["qfcoef"] = vec(qfcoef)
    end
    d["rinner"] = aug.rinner
    isnothing(aug.qijs) || (d["qij"] = [upf_attributes(x, :qij) for x in aug.qijs])
    isnothing(aug.qijls) || (d["qijl"] = [upf_attributes(x, :qijl) for x in aug.qijls])
    return d
end

_from_dict(::Type{UpfChi}, d) = upf_struct(UpfChi, d; chi=upf_get(Vector{Float64}, d, "content"))

upf_wfcs(d, key) =
    [upf_struct(UpfWfc, x; wfc=upf_get(Vector{Float64}, x, "content")) for x in upf_list(d, key)]

_from_dict(::Type{UpfFullWfc}, d) = UpfFullWfc(upf_wfcs(d, "aewfc"), upf_wfcs(d, "pswfc"))

_to_dict(f::UpfFullWfc) = UpfDict("number_of_wfc" => length(f.aewfcs),
                                  "aewfc" => [upf_attributes(w, :wfc) for w in f.aewfcs],
                                  "pswfc" => [upf_attributes(w, :wfc) for w in f.pswfcs])

_from_dict(::Type{UpfSpinOrb}, d) =
    UpfSpinOrb([upf_struct(UpfRelWfc, x) for x in upf_list(d, "relwfc")],
               [upf_struct(UpfRelBeta, x) for x in upf_list(d, "relbeta")])

_to_dict(so::UpfSpinOrb) = UpfDict("relwfc" => upf_attributes.(so.relwfcs),
                                   "relbeta" => upf_attributes.(so.relbetas))

_from_dict(::Type{UpfPaw}, d) =
    upf_struct(UpfPaw, d; aewfcs=upf_wfcs(d, "aewfc"), pswfcs=upf_wfcs(d, "pswfc"))

function _to_dict(paw::UpfPaw)
    d = upf_attributes(paw)
    d["aewfc"] = [upf_attributes(w, :wfc) for w in pop!(d, "aewfcs")]
    d["pswfc"] = [upf_attributes(w, :wfc) for w in pop!(d, "pswfcs")]
    return d
end

function _from_dict(::Type{UpfGipaw}, d)
    orbitals = upf_list(get(d, "gipaw_core_orbitals", nothing), "gipaw_core_orbital")
    return UpfGipaw(upf_get(Int, d, "gipaw_data_format"),
                    [upf_struct(UpfGipawCoreOrbital, x;
                                core_orbital=upf_get(Vector{Float64}, x, "content"))
                     for x in orbitals])
end

function _to_dict(gipaw::UpfGipaw)
    orbitals = [upf_attributes(x, :core_orbital) for x in gipaw.core_orbitals]
    return UpfDict("gipaw_data_format" => gipaw.gipaw_data_format,
                   "gipaw_core_orbitals" => UpfDict("number_of_core_orbitals" => length(orbitals),
                                                    "gipaw_core_orbital" => orbitals))
end
