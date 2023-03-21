@doc raw"""
Analytical Hartwigsen-Goedecker-Hutter pseudopotential.

[C. Hartwigsen, S. Goedecker, and J. Hutter.
*Pys. Rev. B* **58**, 3641 (1998)](https://doi.org/10.1103/PhysRevB.58.3641)
"""
struct HghPsP{T} <: AnalyticalPsP
    "SHA1 Checksum"
    checksum::Vector{UInt8}
    "Atomic charge"
    Zatom::Union{Nothing,Int}
    "Valence charge"
    Zval::Int
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
element(psp::HghPsP) = isnothing(psp.Zatom) ? nothing : PeriodicTable.elements[Int(psp.Zatom)]
has_spin_orbit(::HghPsP)::Bool = false
is_norm_conserving(::HghPsP)::Bool = true
is_ultrasoft(::HghPsP)::Bool = false
is_paw(::HghPsP)::Bool = false
valence_charge(psp::HghPsP) = psp.Zval
atomic_charge(psp::HghPsP) = psp.Zatom
max_angular_momentum(psp::HghPsP)::Int = psp.lmax

has_quantity(::AbstractPsPQuantity, ::HghPsP) = false
has_quantity(::LocalPotential, ::HghPsP) = true
has_quantity(::BetaProjector, ::HghPsP) = true
has_quantity(::BetaCoupling, ::HghPsP) = true

get_quantity(::BetaCoupling, psp::HghPsP) = psp.D
get_quantity(::BetaCoupling, psp::HghPsP, l) = psp.D[l]
get_quantity(::BetaCoupling, psp::HghPsP, l, n) = psp.D[l][n,n]
get_quantity(::BetaCoupling, psp::HghPsP, l, n, m) = psp.D[l][n,m]

n_radials(::BetaProjector, psp::HghPsP, l::Int) = size(psp.D[l], 1)
n_radials(::ChiProjector, ::HghPsP, ::Int) = 0

function cutoff_radius(q::AbstractPsPQuantity, psp::HghPsP; f=nothing, tol=nothing)
    return has_quantity(q, psp) ? nothing : Inf
end
cutoff_radius(q::PsPProjector, psp::HghPsP, l, n; tol=nothing) = cutoff_radius(q, psp)

function psp_quantity_evaluator(::RealSpace, ::LocalPotential, psp::HghPsP)
    return local_potential_real(psp)
end
function psp_quantity_evaluator(::RealSpace, ::BetaProjector, psp::HghPsP, l, n)
    return beta_projector_real(psp, l, n)
end
function psp_quantity_evaluator(::FourierSpace, ::LocalPotential, psp::HghPsP)
    return local_potential_fourier(psp)
end
function psp_quantity_evaluator(::FourierSpace, ::BetaProjector, psp::HghPsP, l, n)
    return beta_projector_fourier(psp, l, n)
end

@doc raw"""
The local potential of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) / (t^2 exp(t^2 / 2))``
where ``x = r_\text{loc} q`` and `Q`
is a polynomial of at most degree 8. This function returns `Q`.
"""
@inline function local_potential_polynomial_fourier(psp::HghPsP, x::T) where {T<:Real}
    rloc::T = psp.rloc
    Zval::T = psp.Zval

    # The polynomial prefactor P(t) (as used inside the { ... } brackets of equation
    # (5) of the HGH98 paper)
    P = (psp.cloc[1]
         + psp.cloc[2] * (3 - x^2)
         + psp.cloc[3] * (15 - 10x^2 + x^4)
         + psp.cloc[4] * (105 - 105x^2 + 21x^4 - x^6))

    return 4T(π) * rloc^2 * (-Zval + sqrt(T(π) / 2) * rloc * x^2 * P)
end

# [GTH98] (6) except they do it with plane waves normalized by 1/sqrt(Ω).
function local_potential_fourier(psp::HghPsP)
    function Vloc(q::T) where {T<:Real}
        x::T = q * psp.rloc
        return local_potential_polynomial_fourier(psp, x) * exp(-x^2 / 2) / x^2
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

# [GTH98] (1)
function local_potential_real(psp::HghPsP)
    function Vloc(r::T)::T where {T<:Real}
        r += iszero(r) ? eps(T) : zero(T)  # quick hack for the division by zero below
        cloc = psp.cloc
        rr = r / psp.rloc
        return -psp.Zval / r * erf(rr / sqrt(T(2))) +
               exp(-rr^2 / 2) * (cloc[1] + cloc[2] * rr^2 + cloc[3] * rr^4 + cloc[4] * rr^6)
    end
    Vloc(R::AbstractVector) = Vloc(norm(R))
    return Vloc
end

@doc raw"""
The nonlocal projectors of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) exp(-t^2 / 2)`` where ``t = r_l q``
and `Q` is a polynomial. This function returns `Q`.
"""
@inline function beta_projector_polynomial_fourier(psp::HghPsP, l::Int, n::Int,
                                              x::T) where {T<:Real}
    @assert 0 <= l <= length(psp.rnl) - 1
    @assert n > 0
    rnl::T = psp.rnl[l]
    common::T = 4T(π)^(5 / T(4)) * sqrt(T(2^(l + 1)) * rnl^3)

    # Note: In the (l == 0 && i == 2) case the HGH paper has an error.
    #       The first 8 in equation (8) should not be under the sqrt-sign
    #       This is the right version (as shown in the GTH paper)
    (l == 0 && n == 1) && return convert(typeof(x), common)
    (l == 0 && n == 2) && return common * 2 / sqrt(T(15)) * (3 - x^2)
    (l == 0 && n == 3) && return common * 4 / 3sqrt(T(105)) * (15 - 10x^2 + x^4)
    #
    (l == 1 && n == 1) && return common * 1 / sqrt(T(3)) * x
    (l == 1 && n == 2) && return common * 2 / sqrt(T(105)) * x * (5 - x^2)
    (l == 1 && n == 3) && return common * 4 / 3sqrt(T(1155)) * x * (35 - 14x^2 + x^4)
    #
    (l == 2 && n == 1) && return common * 1 / sqrt(T(15)) * x^2
    (l == 2 && n == 2) && return common * 2 / 3sqrt(T(105)) * x^2 * (7 - x^2)
    #
    (l == 3 && n == 1) && return common * 1 / sqrt(T(105)) * x^3

    return error("Not implemented for l=$l and i=$n")
end

# [HGH98] (7-15) except they do it with plane waves normalized by 1/sqrt(Ω).
function beta_projector_fourier(psp::HghPsP, l::Int, n::Int)
    function β(q::T) where {T<:Real}
        x::T = q * psp.rnl[l]
        return beta_projector_polynomial_fourier(psp, l, n, x) * exp(-x^2 / T(2))
    end
    β(Q::AbstractVector) = β(norm(Q))
    return β
end

# [HGH98] (3)
function beta_projector_real(psp::HghPsP, l::Int, n::Int)
    function β(r::T) where {T<:Real}
        rnl = T(psp.rnl[l])
        ired = (4n - 1) / T(2)
        return sqrt(T(2)) * r^(l + 2(n - 1)) * exp(-r^2 / 2rnl^2) / rnl^(l + ired) /
               sqrt(gamma(l + ired))
    end
    β(R::AbstractVector) = β(norm(R))
    return β
end

function psp_energy_correction(T::Type{<:Real}, psp::HghPsP)
    # By construction we need to compute the DC component of the difference
    # of the Coulomb potential (-Z/G^2 in Fourier space) and the pseudopotential
    # i.e. -4πZ/(ΔG)^2 -  eval_psp_local_fourier(psp, ΔG) for ΔG → 0. This is:
    cloc_coeffs = [1, 3, 15, 105]
    difference_DC = (T(psp.Zval) * T(psp.rloc)^2 / 2 +
                     sqrt(T(π) / 2) * T(psp.rloc)^3 * T(sum(cloc_coeffs .* psp.cloc)))

    # Multiply by number of electrons and 4π (spherical Hankel prefactor)
    # to get energy per unit cell
    return 4T(π) * difference_DC
end
