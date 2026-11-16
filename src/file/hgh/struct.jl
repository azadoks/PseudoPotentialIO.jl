"""
The header of the HGH formats of ABINIT: `pspcod` 3 (Hartwigsen-Goedecker-Hutter, the
diagonal coupling coefficients) and 10 (Goedecker-Teter-Hutter, all of them).
"""
struct HghAbinitHeader
    "Atomic charge"
    zatom::Float64
    "Pseudo-atomic (valence) charge"
    zion::Float64
    "Generation date `ddmmyy`"
    pspdat::Int
    "PSP format: 3 or 10"
    pspcod::Int
    "ABINIT code for exchange-correlation if positive, LibXC code if negative"
    pspxc::Int
    "Angular momentum channel of the local potential (unused)"
    lloc::Int
    "Unused"
    mmax::Int
    "Unused"
    r2well::Float64
    HghAbinitHeader(fields...) = validated(new(fields...))
end

function validate(h::HghAbinitHeader)
    @check "HGH" h.pspcod in (3, 10) "pspcod = $(h.pspcod) is not an HGH format (3 or 10)"
    @check "HGH" 0 < h.zion <= h.zatom "zion = $(h.zion), zatom = $(h.zatom)"
    @check "HGH" psp8_known_pspxc(h.pspxc) "pspxc = $(h.pspxc) is not a known functional"
end

"""
Hartwigsen-Goedecker-Hutter pseudopotential file contents: the CP2K format, or
the HGH formats of ABINIT (`pspcod` 3 and 10, with `abinit`).
"""
struct HghFile <: PsPFile
    "Identifier"
    identifier::String
    "Functional: from `pspxc` (ABINIT formats), or guessed from the source file path"
    functional::String
    "Description"
    title::String
    "Number of valence electrons of each angular momentum (CP2K format; empty for the ABINIT formats, whose total is `abinit.zion`)"
    zion::Vector{Int}
    "Cutoff radius for the local part of the pseudopotential"
    rloc::Float64
    "Number of coefficients defining the local part of the pseudopotential"
    nloc::Int
    "Coefficients of the local part of the pseudopotential"
    cloc::Vector{Float64}
    "Maximum angular momentum"
    lmax::Int
    "Non-local projector cutoff radius for each angular momentum"
    rp::Vector{Float64}
    "Kleinman-Bylander energies"
    h::Vector{Matrix{Float64}}
    "Spin-orbit coupling matrices, for each angular momentum (0×0 for l = 0, and in the CP2K format, which has none)"
    k::Vector{Matrix{Float64}}
    "Header of the ABINIT formats; `nothing` for the CP2K format"
    abinit::Union{Nothing,HghAbinitHeader}
    HghFile(fields...) = validated(new(fields...))
end

#= The coupling coefficients of ABINIT's `pspcod = 3` (HGH 1998, Eqs. (6)-(8)): only the
diagonal is in the file; the off-diagonal follow from it, for l = 0, 1, 2, as the factors
(h12 / h22, h13 / h33, h23 / h33). l = 3 has one projector. =#
const HGH_OFFDIAGONAL = ((-sqrt(3 / 5) / 2, sqrt(5 / 21) / 2, -sqrt(100 / 63) / 2),
                         (-sqrt(5 / 7) / 2, sqrt(35 / 11) / 6, -14 / sqrt(11) / 6),
                         (-sqrt(7 / 9) / 2, sqrt(63 / 143) / 2, -18 / sqrt(143) / 2))

hgh_max_projectors(l) = l < 3 ? 3 : 1

# The coupling matrix of angular momentum l from its diagonal (`pspcod = 3`)
function hgh_coupling_matrix(l::Integer, diagonal::AbstractVector)
    n = length(diagonal)
    n <= hgh_max_projectors(l) ||
        error("$n projectors for l = $l, but the HGH format has at most $(hgh_max_projectors(l))")
    h = zeros(n, n)
    for i in 1:n
        h[i, i] = diagonal[i]
    end
    if n >= 2
        h[1, 2] = h[2, 1] = HGH_OFFDIAGONAL[l + 1][1] * diagonal[2]
    end
    if n >= 3
        h[1, 3] = h[3, 1] = HGH_OFFDIAGONAL[l + 1][2] * diagonal[3]
        h[2, 3] = h[3, 2] = HGH_OFFDIAGONAL[l + 1][3] * diagonal[3]
    end
    return h
end

hgh_diagonal(h::AbstractMatrix) = [h[i, i] for i in axes(h, 1)]

function validate(file::HghFile)
    header = file.abinit
    if isnothing(header)
        @check "HGH" all(>=(0), file.zion) && sum(file.zion) > 0 "zion = $(file.zion)"
        @check "HGH" all(isempty, file.k) "spin-orbit coupling (k), which the CP2K format does not have"
    else
        @check "HGH" isempty(file.zion) "zion = $(file.zion) for an ABINIT format (its zion is in the header)"
        @check "HGH" 0 <= file.lmax <= 3 "lmax = $(file.lmax) is not in 0:3"
        if header.pspcod == 3
            @check "HGH" file.nloc == 4 "pspcod = 3 has 4 local coefficients, not $(file.nloc)"
        end
    end
    @check "HGH" file.rloc > 0 "rloc = $(file.rloc)"
    @check "HGH" file.nloc == length(file.cloc) "nloc = $(file.nloc), but $(length(file.cloc)) coefficients"
    @check "HGH" 0 <= file.nloc <= 4 "nloc = $(file.nloc) is not in 0:4"
    @check "HGH" length(file.rp) == length(file.h) == length(file.k) == file.lmax + 1 "lmax = $(file.lmax), but $(length(file.rp)) radii, $(length(file.h)) coupling matrices and $(length(file.k)) spin-orbit ones"
    for (l, (rp, h, k)) in zip(0:file.lmax, zip(file.rp, file.h, file.k))
        @check "HGH" rp > 0 || (isempty(h) && isempty(k)) "rp = $rp for l = $l"
        for (name, m) in (("coupling", h), ("spin-orbit coupling", k))
            @check "HGH" size(m, 1) == size(m, 2) <= 3 && m == permutedims(m) "the $name matrix of l = $l is not symmetric of size ≤ 3"
        end
        @check "HGH" l > 0 || isempty(k) "spin-orbit coupling for l = 0"
        isnothing(header) && continue
        if header.pspcod == 3
            for (name, m) in (("coupling", h), ("spin-orbit coupling", k))
                @check "HGH" size(m, 1) <= hgh_max_projectors(l) && m ≈ hgh_coupling_matrix(l, hgh_diagonal(m)) "the $name matrix of l = $l does not follow the HGH relations of pspcod = 3"
            end
        elseif l > 0
            @check "HGH" size(k) == size(h) "pspcod = 10: for l = $l, $(size(h, 1)) projectors but $(size(k, 1)) spin-orbit ones"
        end
    end
end

"""
Read an HGH file (its text) into nested dictionaries: of the CP2K format
(`src/file/hgh/file_dict_cp2k.jl`) or of an ABINIT format (`src/file/hgh/file_dict_abinit.jl`).
"""
hgh_to_dict(text::AbstractString) =
    hgh_is_abinit(text) ? hgh_abinit_to_dict(text) : hgh_cp2k_to_dict(text)

"""
Write nested dictionaries (see `hgh_to_dict`) as an HGH file to `io`, in the format they are
of (ABINIT's if they have a "pspcod").
"""
dict_to_hgh(io::IO, d::AbstractDict) =
    haskey(d, "pspcod") ? dict_to_hgh_abinit(io, d) : dict_to_hgh_cp2k(io, d)

HghFile(io::IO; identifier="", validation=VALIDATION[]) =
    with_validation(validation; source=identifier) do
        return _from_dict(HghFile, hgh_to_dict(read(io, String)); identifier)
    end

function HghFile(path::AbstractString; identifier="", validation=VALIDATION[])
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return HghFile(io; identifier, validation)
    end
end

identifier(psp::HghFile)::String = psp.identifier
format(::HghFile)::String = "HGH"
function element(psp::HghFile)::PeriodicTable.Element
    isnothing(psp.abinit) || return PeriodicTable.elements[Int(psp.abinit.zatom)]
    # Guess the element from the title line, if possible
    title = split(psp.title)
    isempty(title) && return PeriodicTable.Element()
    symbol = title[1]
    return haskey(PeriodicTable.elements, Symbol(symbol)) ? PeriodicTable.elements[Symbol(symbol)] : PeriodicTable.Element()
end
functional(psp::HghFile)::Vector{Symbol} = map(name -> Symbol(lowercase(name)), split(psp.functional))
function valence_charge(psp::HghFile)::Float64
    return isnothing(psp.abinit) ? sum(psp.zion; init=0) : psp.abinit.zion
end
is_norm_conserving(::HghFile)::Bool = true
is_ultrasoft(::HghFile)::Bool = false
is_paw(::HghFile)::Bool = false
has_spin_orbit(psp::HghFile)::Bool = any(k -> any(!iszero, k), psp.k)
has_model_core_charge_density(::HghFile)::Bool = false

save_psp(io::IO, file::HghFile) = dict_to_hgh(io, _to_dict(file))
