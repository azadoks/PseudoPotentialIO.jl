#= PSP8 dictionaries (`file_dict.jl`) <-> `Psp8File`: `_from_dict` builds the struct, with
the sizes and an entry (empty if there is no block) for every l in `projectors` and `ekb`,
and for every l ≥ 1 in their spin-orbit counterparts (l = 0 is empty); `_to_dict` is its
inverse. =#

# The projectors and ekb of the blocks `blocks`, by l in 0:lmax
function psp8_by_l(blocks, lmax)
    projectors = [Vector{Float64}[] for _ in 0:lmax]
    ekb = [Float64[] for _ in 0:lmax]
    for block in blocks
        projectors[block["l"] + 1] = block["projectors"]
        ekb[block["l"] + 1] = block["ekb"]
    end
    return projectors, ekb
end

function _from_dict(::Type{Psp8File}, d::AbstractDict; identifier="")
    h = d["header"]
    lmax = h["lmax"]
    projectors, ekb = psp8_by_l(d["projectors"], lmax)
    projectors_so, ekb_so = haskey(d, "projectors_so") ?
                            psp8_by_l(d["projectors_so"], lmax) : (nothing, nothing)
    nprojso = isnothing(projectors_so) ? nothing : length.(projectors_so)
    header = Psp8Header(h["title"], h["zatom"], h["zion"], h["pspd"], h["pspcod"],
                        h["pspxc"], lmax, h["lloc"], length(d["rgrid"]), h["r2well"],
                        h["rchrg"], h["fchrg"], h["qchrg"], length.(projectors),
                        h["extension_switch"], nprojso)
    nlcc = get(d, "nlcc", Dict())
    valence = get(d, "valence", Dict())
    return Psp8File(identifier, header, d["rgrid"], d["v_local"], projectors, ekb,
                    projectors_so, ekb_so,
                    (get(nlcc, k, nothing) for k in PSP8_NLCC_COLUMNS)...,
                    (get(valence, k, nothing) for k in PSP8_VALENCE_COLUMNS)...)
end

# the blocks of the l that have projectors
psp8_blocks(projectors, ekb) =
    [Psp8Dict("l" => l - 1, "ekb" => ekb[l], "projectors" => projectors[l])
     for l in eachindex(projectors) if !isempty(projectors[l])]

function _to_dict(file::Psp8File)
    h = file.header
    header = Psp8Dict(string(f) => getfield(h, f) for f in fieldnames(Psp8Header)
                      if !(f in (:mmax, :nproj, :nprojso)))
    d = Psp8Dict("header" => header, "rgrid" => file.rgrid,
                 "projectors" => psp8_blocks(file.projectors, file.ekb),
                 "v_local" => file.v_local)
    isnothing(file.projectors_so) ||
        (d["projectors_so"] = psp8_blocks(file.projectors_so, file.ekb_so))
    isnothing(file.rhoc) ||
        (d["nlcc"] = Psp8Dict(k => getfield(file, Symbol(k)) for k in PSP8_NLCC_COLUMNS))
    if !isnothing(file.rhov)
        columns = [getfield(file, Symbol(k)) for k in PSP8_VALENCE_COLUMNS]
        # a column after a missing one is written after zeros
        !isnothing(columns[3]) && isnothing(columns[2]) && (columns[2] = zero(file.rhov))
        d["valence"] = Psp8Dict(k => c for (k, c) in zip(PSP8_VALENCE_COLUMNS, columns)
                                if !isnothing(c))
    end
    return d
end
