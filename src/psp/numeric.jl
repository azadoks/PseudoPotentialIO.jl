@doc raw"""
Abstract type representing numeric pseudopotentials.

All quantities must be in Hartree atomic units without prefactors like `r` or `4π`.

- Lengths in Bohr radii (a₀)
- Energies in Hartree (Ha / Eₕ)
- Electric charge in electron charges (e = 1)
- Mass in electron masses (mₑ)
- Action in reduced Plank constants (ħ = 1)

Vectors indexed by angular momentum should be `OffsetVector`s with indices starting
at zero so that angular momentum `l`, which naturally starts at 0, can be used for
both computation and indexing.

Required fields:

```julia
# Atomic total charge in units of electron charge
Ztot::Number
# Pseudo-atomic valence charge in units of electron charge
Zval::Number
# Maximum angular momentum
lmax::Integer
# Radial mesh in units of Bohr
r::AbstractVector{Real}
# Radial mesh spacing in units of Bohr
dr::Union{Real, AbstractVector{Real}}
# Local potential on the radial mesh in units of Hartree
Vloc::AbstractVector{Real}
## The units of `D` and `β` should be such that `⟨ βˡₙ | Dˡₙₘ | βˡₘ ⟩` gives Hartree
# Nonlocal projector coupling constants D[l][n,m]
D::OffsetVector{AbstractMatrix{Real}}
# Nonlocal projectors on the radial mesh β[l][n]
β::OffsetVector{AbstractVector{AbstractVector{Real}}}
# Cutoff indices for nonlocal projectors
β_ircut::OffsetVector{AbstractVector{Int}}

## "Optional" fields (must still exist, but could be Union{Nothing})
# Model core charge density (non-linear core correction) on the radial mesh
ρcore::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic valence charge density on the radial mesh
ρval::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic orbitals on the radial mesh ϕ̃[l][n]
ϕ̃::Union{Nothing,OffsetVector{AbstractVector{AbstractVector{Real}}}}
# Cutoff indices for pseudo-atomic orbitals
ϕ̃_ircut::Union{Nothing,OffsetVector{AbstractVector{Int}}}
```
"""
abstract type NumericPsP{T} <: AbstractPsP end

function element(psp::NumericPsP)::String
    return PeriodicTable.elements[Int(psp.Ztot)].symbol
end
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
n_projectors(psp::NumericPsP, l::Int)::Int = length(psp.β[l])
n_pseudo_orbitals(psp::NumericPsP, l::Int)::Int = isnothing(psp.ϕ̃) ? 0 : length(psp.ϕ̃[l])
valence_charge(psp::NumericPsP{T})::T = psp.Zval
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation
has_nlcc(psp::NumericPsP)::Bool = !isnothing(psp.ρcore)

function projector_coupling(psp::NumericPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.D[l]
end

#TODO test the *_real functions
#TODO test `nothing` cases

function projector_real(psp::NumericPsP, l::Int, n::Int, r::T)::T where {T<:Real}
    return interpolate((psp.r,), psp.β[l][n], (Gridded(Linear()),))(r)
end

function pseudo_orbital_real(psp::NumericPsP, l::Int, n::Int,
                             r::T)::Union{Nothing,T} where {T<:Real}
    isnothing(psp.ϕ̃) && return nothing
    return interpolate((psp.r,), psp.ϕ̃[l][n], (Gridded(Linear()),))(r)
end

function valence_charge_density_real(psp::NumericPsP,
                                     r::T)::Union{Nothing,T} where {T<:Real}
    isnothing(psp.ρval) && return nothing
    return interpolate((psp.r,), psp.ρval, (Gridded(Linear()),))(r)
end

function core_charge_density_real(psp::NumericPsP, r::T)::Union{Nothing,T} where {T<:Real}
    !has_nlcc(psp) && return nothing
    return interpolate((psp.r,), psp.ρcore, (Gridded(Linear()),))(r)
end

function local_potential_fourier(psp::NumericPsP, q::T)::T where {T<:Real}
    f = @. psp.r * fast_sphericalbesselj0(q * psp.r) * (psp.r * psp.Vloc + psp.Zval)
    return 4π * (simpson(f, psp.dr) - psp.Zval / q^2)
end

function projector_fourier(psp::NumericPsP, l::Int, n::Int, q::T)::T where {T<:Real}
    r = @view psp.r[1:psp.β_ircut[l][n]]
    f = @. r^2 * fast_sphericalbesselj(l, q * r) * psp.β[l][n]
    return 4π * trapezoid(f, psp.dr)
end

function pseudo_orbital_fourier(psp::NumericPsP, l::Int, n::Int,
                                q::T)::Union{Nothing,T} where {T<:Real}
    isnothing(psp.ϕ̃) && return nothing
    r = @view psp.r[1:psp.ϕ̃_ircut[l][n]]
    f = @. r^2 * fast_sphericalbesselj(l, q * r) * psp.ϕ̃[l][n]
    return 4π * trapezoid(f, psp.dr)
end

function valence_charge_density_fourier(psp::NumericPsP,
                                        q::T)::Union{Nothing,T} where {T<:Real}
    isnothing(psp.ρval) && return nothing
    f = @. psp.r^2 * fast_sphericalbesselj0(q * psp.r) * psp.ρval
    return 4π * trapezoid(f, psp.dr)
end

function core_charge_density_fourier(psp::NumericPsP,
                                     q::T)::Union{Nothing,T} where {T<:Real}
    isnothing(psp.ρcore) && return nothing
    f = @. psp.r^2 * fast_sphericalbesselj0(q * psp.r) * psp.ρcore
    return 4π * trapezoid(f, psp.dr)
end

function pseudo_energy_correction(psp::NumericPsP{T})::T where {T<:Real}
    f = @. psp.r * (psp.r * psp.Vloc + psp.Zval)
    return 4π * trapezoid(f, psp.dr)
end