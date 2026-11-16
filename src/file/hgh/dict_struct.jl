#= HGH dictionaries (`file_dict.jl`) <-> `HghFile`: `_from_dict` builds the struct, with the
sizes and the symmetric coupling matrices; `_to_dict` is its inverse. =#

function _from_dict(::Type{HghFile}, d::AbstractDict; identifier="")
    haskey(d, "pspcod") && return hgh_abinit_from_dict(d; identifier)
    # the file does not say: guess from its name
    if occursin("lda", identifier)
        functional = "LDA_C LDA_X_PZ"
    elseif occursin("pbe", identifier)
        functional = "GGA_C_PBE GGA_X_PBE"
    else
        functional = ""
    end
    nonlocal = d["nonlocal"]
    return HghFile(identifier, functional, d["title"], d["zion"], d["rloc"],
                   length(d["cloc"]), d["cloc"], length(nonlocal) - 1,
                   [channel["rp"] for channel in nonlocal],
                   [hgh_matrix(channel["h"]) for channel in nonlocal],
                   [zeros(0, 0) for _ in nonlocal], nothing)
end

function hgh_abinit_from_dict(d::AbstractDict; identifier="")
    pspcod, nonlocal = d["pspcod"], d["nonlocal"]
    lmax = length(nonlocal) - 1
    @check "HGH" d["lmax"] == lmax "lmax = $(d["lmax"]), but $(lmax + 1) angular momenta"
    k = map(0:lmax, nonlocal) do l, channel
        l == 0 && return zeros(0, 0)
        return pspcod == 3 ? hgh_coupling_matrix(l, hgh_trim(channel["k"])) : hgh_matrix(channel["k"])
    end
    h = map(0:lmax, nonlocal) do l, channel
        return pspcod == 3 ? hgh_coupling_matrix(l, hgh_trim(channel["h"])) : hgh_matrix(channel["h"])
    end
    header = HghAbinitHeader(d["zatom"], d["zion"], d["pspdat"], pspcod, d["pspxc"], d["lloc"],
                             d["mmax"], d["r2well"])
    return HghFile(identifier, abinit_libxc_string(d["pspxc"]), d["title"], Int[], d["rloc"],
                   length(d["cloc"]), d["cloc"], lmax, [channel["rp"] for channel in nonlocal],
                   h, k, header)
end

# The symmetric matrix of the rows of its upper triangle
function hgh_matrix(rows)
    n = length(rows)
    return [rows[min(i, j)][abs(j - i) + 1] for i in 1:n, j in 1:n]
end

# The rows of the upper triangle of a matrix
hgh_rows(m::AbstractMatrix) = [m[i, i:end] for i in axes(m, 1)]

# `pspcod = 3`: the projectors are those up to the last nonzero diagonal coefficient
hgh_trim(diagonal) = diagonal[1:something(findlast(!iszero, diagonal), 0)]

function _to_dict(file::HghFile)
    isnothing(file.abinit) || return hgh_abinit_to_dict(file)
    nonlocal = [HghDict("rp" => rp, "h" => hgh_rows(h)) for (rp, h) in zip(file.rp, file.h)]
    return HghDict("title" => file.title, "zion" => file.zion, "rloc" => file.rloc,
                   "cloc" => file.cloc, "nonlocal" => nonlocal)
end

function hgh_abinit_to_dict(file::HghFile)
    header = file.abinit
    pspcod = header.pspcod
    # pspcod = 3: the diagonal, with zeros up to the projectors of the format
    padded(l, m) = [hgh_diagonal(m); zeros(hgh_max_projectors(l) - size(m, 1))]
    nonlocal = map(0:file.lmax, file.rp, file.h, file.k) do l, rp, h, k
        channel = HghDict("rp" => rp, "h" => pspcod == 3 ? padded(l, h) : hgh_rows(h))
        l > 0 && (channel["k"] = pspcod == 3 ? padded(l, k) : hgh_rows(k))
        return channel
    end
    return HghDict("title" => file.title, "zatom" => header.zatom, "zion" => header.zion,
                   "pspdat" => header.pspdat, "pspcod" => pspcod, "pspxc" => header.pspxc,
                   "lmax" => file.lmax, "lloc" => header.lloc, "mmax" => header.mmax,
                   "r2well" => header.r2well, "rloc" => file.rloc, "cloc" => file.cloc,
                   "nonlocal" => nonlocal)
end
