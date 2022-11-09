@doc raw"""
Abstract type representing numeric pseudopotentials.

All quantities must be in atomic units without prefactors like `r` or `4π`.
Vectors indexed by angular momentum should be `OffsetVector`s with indices starting
at zero.

Required fields:

```julia
# Atomic total charge
Ztot::Number
# Pseudo-atomic valence charge
Zval::Number
# Maximum angular momentum
lmax::Integer
# Radial mesh
r::AbstractVector{Real}
# Radial mesh spacing
dr::Union{Real, AbstractVector{Real}}
# Local potential on the radial mesh
Vloc::AbstractVector{Real}
# Nonlocal projector coupling constants D[l][n,m]
D::OffsetVector{AbstractMatrix{Real}}
# Nonlocal projectors on the radial mesh β[l][n]
β::OffsetVector{AbstractVector{AbstractVector{Real}}}
# Cutoff indices for nonlocal projectors
β_ircut::OffsetVector{AbstractVector{Int}}
```

Optional fields:

```julia
# Model core charge density (non-linear core correction) on the radial mesh
ρcore::AbstractVector{Real}
# Pseudo-atomic valence charge density on the radial mesh
ρval::AbstractVector{Real}
# Pseudo-atomic orbitals on the radial mesh ϕ̃[l][n]
ϕ̃::OffsetVector{AbstractVector{AbstractVector{Real}}}
# Cutoff indices for pseudo-atomic orbitals
ϕ̃_ircut::OffsetVector{AbstractVector{Int}}
```
"""
abstract type NumericPsP{T} <: AbstractPsP end

function element(psp::NumericPsP)::PeriodicTable.Element
    return PeriodicTable.elements[Int(psp.Ztot)]
end

relativistic_treatment(::NumericPsP)::Symbol = :scalar
has_nlcc(psp::NumericPsP)::Bool = !isnothing(psp.ρcore)
n_projectors(psp::NumericPsP)::Int = sum(length.(psp.β))
n_pseudo_orbitals(psp::NumericPsP)::Int = isnothing(psp.ϕ̃) ? 0 : sum(length.(psp.β))

function max_angular_momentum(psp::NumericPsP)::Int
    return psp.lmax
end

function valence_charge(psp::NumericPsP{T})::T where {T<:Real}
    return psp.Zval
end

function projector_coupling(psp::NumericPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.D[l]
end

function projector_coupling(psp::NumericPsP{T}, l::Int, n::Int)::T where {T<:Real}
    return psp.D[l][n, n]
end

function projector_coupling(psp::NumericPsP{T}, l::Int, ni::Int, nj::Int)::T where {T<:Real}
    return psp.D[l][ni, nj]
end

function local_potential_fourier(psp::NumericPsP, q::T)::T where {T<:Real}
    f = @. psp.r * fast_sphericalbesselj0(q * psp.r) * (psp.r * psp.Vloc + psp.Zval)
    return 4π * (simpson(f, psp.dr) - psp.Zval / q^2)
end

function pseudo_energy_correction(psp::NumericPsP{T})::T where {T<:Real}
    f = @. psp.r * (psp.r * psp.Vloc + psp.Zval)
    4π * trapezoid(f, psp.dr)
end

function projector_fourier(psp::NumericPsP, l::Int, n::Int, q::T)::T where {T<:Real}
    r = @view psp.r[1:psp.β_ircut[l][n]]
    f = @. r^2 * fast_sphericalbesselj(l, q * r) * psp.β[l][n]
    return 4π * trapezoid(f, psp.dr)
end

function core_charge_density_fourier(psp::NumericPsP, q::T)::T where {T<:Real}
    f = @. psp.r^2 * fast_sphericalbesselj0(q * psp.r) * psp.ρcore
    return 4π * trapezoid(f, psp.dr)
end

function valence_charge_density_fourier(psp::NumericPsP, q::T)::T where {T<:Real}
    f = @. psp.r^2 * fast_sphericalbesselj0(q * psp.r) * psp.ρval
    return 4π * trapezoid(f, psp.dr)
end

function pseudo_orbital_fourier(psp::NumericPsP, l::Int, n::Int, q::T)::T where {T<:Real}
    r = @view psp.r[1:psp.ϕ̃_ircut[l][n]]
    f = @. r^2 * fast_sphericalbesselj(l, q * r) * psp.ϕ̃[l][n]
    return 4π * trapezoid(f, psp.dr)
end
