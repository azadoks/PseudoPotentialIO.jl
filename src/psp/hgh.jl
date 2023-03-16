@doc raw"""
Analytical Hartwigsen-Goedecker-Hutter pseudopotential.

[C. Hartwigsen, S. Goedecker, and J. Hutter.
*Pys. Rev. B* **58**, 3641 (1998)](https://doi.org/10.1103/PhysRevB.58.3641)
"""
struct HghPsP{T} <: AnalyticalPsP
    "SHA1 Checksum"
    checksum::Vector{UInt8}
    "Atomic charge"
    Zatom::Union{Nothing,T}
    "Valence charge"
    Zval::T
    "Maximum angular momentum"
    lmax::Int
    "Radial cutoff for the local part of the pseudopotential"
    rloc::T
    "Polynomial coefficience of the local part of the pseudopotential"
    cloc::Vector{T}
    "Radial cutoffs for the nonlocal projectors `rnl[l]`"
    rnl::OffsetVector{T,Vector{T}}
    "Nonlocal projector coupling coefficients `D[l][n,m]`"
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
end

function HghPsP(file::HghFile)
    el = element(file)
    Zatom = isnothing(el) ? nothing : PeriodicTable.elements[Symbol(el)].number
    cloc = file.cloc
    length(cloc) <= 4 || error("length(cloc) > 4 not supported.")
    if length(cloc) < 4
        n_extra = 4 - length(cloc)
        cloc = [cloc; zeros(Float64, n_extra)]
    end

    rnl = OffsetVector(file.rp, 0:(file.lmax))
    D = OffsetVector(file.h, 0:(file.lmax))

    return HghPsP{Float64}(file.checksum, Zatom, sum(Float64.(file.zion)), file.lmax,
                           file.rloc, cloc, rnl, D)
end

identifier(psp::HghPsP)::String = bytes2hex(psp.checksum)
element(psp::HghPsP)::String = isnothing(psp.Zatom) ? "unknown" :
                               PeriodicTable.elements[Int(psp.Zatom)].symbol
has_spin_orbit(::HghPsP)::Bool = false
has_core_density(::HghPsP)::Bool = false
has_valence_density(::HghPsP)::Bool = false
has_pseudo_orbitals(::HghPsP)::Bool = false
is_norm_conserving(::HghPsP)::Bool = true
is_ultrasoft(::HghPsP)::Bool = false
is_paw(::HghPsP)::Bool = false
valence_charge(psp::HghPsP) = psp.Zval
atomic_charge(psp::HghPsP) = psp.Zatom
max_angular_momentum(psp::HghPsP)::Int = psp.lmax
n_projector_radials(psp::HghPsP, l::Int)::Int = size(psp.D[l], 1)
n_pseudo_orbital_radials(::HghPsP, ::Int)::Int = 0

local_potential_cutoff_radius(::HghPsP; tol=nothing) = Inf
projector_cutoff_radius(::HghPsP, ::Int, ::Int; tol=nothing) = Inf
pseudo_orbital_cutoff_radius(::HghPsP, ::Int, ::Int; tol=nothing) = Inf
valence_charge_density_cutoff_radius(::HghPsP; tol=nothing) = Inf
core_charge_density_cutoff_radius(::HghPsP; tol=nothing) = Inf

projector_coupling(psp::HghPsP, l::Int) = psp.D[l]

@doc raw"""
The local potential of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) / (t^2 exp(t^2 / 2))``
where ``x = r_\text{loc} q`` and `Q`
is a polynomial of at most degree 8. This function returns `Q`.
"""
@inline function local_potential_polynomial_fourier(psp::HghPsP, x::T)::T where {T<:Real}
    rloc = psp.rloc
    Zval = psp.Zval

    # The polynomial prefactor P(t) (as used inside the { ... } brackets of equation
    # (5) of the HGH98 paper)
    P = (psp.cloc[1]
         + psp.cloc[2] * (3 - x^2)
         + psp.cloc[3] * (15 - 10x^2 + x^4)
         + psp.cloc[4] * (105 - 105x^2 + 21x^4 - x^6))

    return 4π * rloc^2 * (-Zval + sqrt(π / 2) * rloc * x^2 * P)
end

# [GTH98] (6) except they do it with plane waves normalized by 1/sqrt(Ω).
function local_potential_fourier(psp::HghPsP)
    function Vloc(q)
        x = q * psp.rloc
        return local_potential_polynomial_fourier(psp, x) * exp(-x^2 / 2) / x^2
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

# [GTH98] (1)
function local_potential_real(psp::HghPsP)
    function Vloc(r::T) where {T}
        r += iszero(r) ? eps(T) : zero(T)  # quick hack for the division by zero below
        cloc = psp.cloc
        rr = r / psp.rloc
        return convert(T,
                       -psp.Zval / r * erf(rr / sqrt(T(2))) +
                       exp(-rr^2 / 2) *
                       (cloc[1] + cloc[2] * rr^2 + cloc[3] * rr^4 + cloc[4] * rr^6))
    end
    Vloc(R::AbstractVector) = Vloc(norm(R))
    return Vloc
end

@doc raw"""
The nonlocal projectors of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) exp(-t^2 / 2)`` where ``t = r_l q``
and `Q` is a polynomial. This function returns `Q`.
"""
@inline function projector_polynomial_fourier(psp::HghPsP, l::Int, n::Int,
                                              x::T)::T where {T<:Real}
    @assert 0 <= l <= length(psp.rnl) - 1
    @assert n > 0
    rnl = psp.rnl[l]
    common = 4π^(5 / 4) * sqrt(2^(l + 1) * rnl^3)

    # Note: In the (l == 0 && i == 2) case the HGH paper has an error.
    #       The first 8 in equation (8) should not be under the sqrt-sign
    #       This is the right version (as shown in the GTH paper)
    (l == 0 && n == 1) && return convert(typeof(x), common)
    (l == 0 && n == 2) && return common * 2 / sqrt(15) * (3 - x^2)
    (l == 0 && n == 3) && return common * 4 / 3sqrt(105) * (15 - 10x^2 + x^4)
    #
    (l == 1 && n == 1) && return common * 1 / sqrt(3) * x
    (l == 1 && n == 2) && return common * 2 / sqrt(105) * x * (5 - x^2)
    (l == 1 && n == 3) && return common * 4 / 3sqrt(1155) * x * (35 - 14x^2 + x^4)
    #
    (l == 2 && n == 1) && return common * 1 / sqrt(15) * x^2
    (l == 2 && n == 2) && return common * 2 / 3sqrt(105) * x^2 * (7 - x^2)
    #
    (l == 3 && n == 1) && return common * 1 / sqrt(105) * x^3

    return error("Not implemented for l=$l and i=$n")
end

# [HGH98] (7-15) except they do it with plane waves normalized by 1/sqrt(Ω).
function projector_fourier(psp::HghPsP, l::Int, n::Int)
    function β(q::T) where {T}
        x::T = q * psp.rnl[l]
        return projector_polynomial_fourier(psp, l, n, x) * exp(-x^2 / 2)
    end
    β(Q::AbstractVector) = β(norm(Q))
    return β
end

# [HGH98] (3)
function projector_real(psp::HghPsP, l::Int, n::Int)
    function β(r::T) where {T}
        rnl = T(psp.rnl[l])
        ired = (4n - 1) / T(2)
        return sqrt(T(2)) * r^(l + 2(n - 1)) * exp(-r^2 / 2rnl^2) / rnl^(l + ired) /
               sqrt(gamma(l + ired))
    end
    β(R::AbstractVector) = β(norm(R))
    return β
end

function pseudo_energy_correction(T::Type, psp::HghPsP)
    # By construction we need to compute the DC component of the difference
    # of the Coulomb potential (-Z/G^2 in Fourier space) and the pseudopotential
    # i.e. -4πZ/(ΔG)^2 -  eval_psp_local_fourier(psp, ΔG) for ΔG → 0. This is:
    cloc_coeffs = [1, 3, 15, 105]
    difference_DC = (psp.Zval * psp.rloc^2 / 2 +
                     sqrt(π / 2) * psp.rloc^3 * sum(cloc_coeffs .* psp.cloc))

    # Multiply by number of electrons and 4π (spherical Hankel prefactor)
    # to get energy per unit cell
    return T(4π * difference_DC)
end
