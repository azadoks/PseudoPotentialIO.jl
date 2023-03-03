@doc raw"""
Abstract type representing numeric pseudopotentials.

All quantities must be in Hartree atomic units.

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
# Nonlocal projectors on the radial mesh, multiplied by the mesh squared r^2 β[l][n]
β::OffsetVector{AbstractVector{AbstractVector{Real}}}

## "Optional" fields (must still exist, but could be Union{Nothing})
# Model core charge density (non-linear core correction) on the radial mesh, multiplied by
# the mesh squared r^2 ρcore
ρcore::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic valence charge density on the radial mesh, multiplied by the mesh squared
# r^2 ρval
ρval::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic orbitals on the radial mesh, multiplied by the mesh squared r^2 ϕ̃[l][n]
ϕ̃::Union{Nothing,OffsetVector{AbstractVector{AbstractVector{Real}}}}
```
"""
abstract type NumericPsP{T} <: AbstractPsP end

function element(psp::NumericPsP)::String
    return PeriodicTable.elements[Int(psp.Ztot)].symbol
end
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
n_projector_radials(psp::NumericPsP, l::Int)::Int = length(psp.β[l])
n_pseudo_orbital_radials(psp::NumericPsP, l::Int)::Int = isnothing(psp.ϕ̃) ? 0 : length(psp.ϕ̃[l])
valence_charge(psp::NumericPsP) = psp.Zval
atomic_charge(psp::NumericPsP) = psp.Ztot
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation
has_nlcc(psp::NumericPsP)::Bool = !isnothing(psp.ρcore)
has_ρval(psp::NumericPsP)::Bool = !isnothing(psp.ρval)
has_ϕ̃(psp::NumericPsP)::Bool = !isnothing(psp.ϕ̃)

function projector_coupling(psp::NumericPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.D[l]
end

function local_potential_real(psp::NumericPsP)
    itp = build_interpolator(psp.Vloc, psp.r)
    Vloc(r) = itp(r)
    Vloc(R::AbstractVector) = Vloc(norm(R))
    return Vloc
end

function projector_real(psp::NumericPsP, l::Int, n::Int)
    isnothing(psp.β[l][n]) && return _ -> nothing
    itp = build_interpolator(psp.β[l][n], psp.r)
    β(r) = itp(r)
    β(R::AbstractVector) = β(norm(R))
    return β
end

function pseudo_orbital_real(psp::NumericPsP, l::Int, n::Int)
    isnothing(psp.ϕ̃) && return _ -> nothing
    itp = build_interpolator(psp.ϕ̃[l][n], psp.r)
    ϕ̃(r) = itp(r)
    ϕ̃(R::AbstractVector) = ϕ̃(norm(R))
    return ϕ̃
end

function valence_charge_density_real(psp::NumericPsP)
    isnothing(psp.ρval) && return _ -> nothing
    itp = build_interpolator(psp.ρval, psp.r)
    ρval(r) = itp(r)
    ρval(R::AbstractVector) = ρval(norm(R))
    return ρval
end

function core_charge_density_real(psp::NumericPsP)
    isnothing(psp.ρcore) && return _ -> nothing
    itp = build_interpolator(psp.ρcore, psp.r)
    ρcore(r) = itp(r)
    ρval(R::AbstractVector) = ρcore(norm(R))
end

@inbounds function local_potential_fourier(psp::NumericPsP; tol=nothing)
    i_start = firstindex(psp.Vloc)
    i_stop = find_truncation_index(psp.Vloc, tol)
    function Vloc(q)
        integrand(i::Int) = psp.r[i] * fast_sphericalbesselj0(q * psp.r[i]) *
            (psp.r[i] * psp.Vloc[i] + psp.Zval)
        integral = dotprod(integrand, i_start, i_stop, psp.dr)
        4π * (integral - psp.Zval / q^2)
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

@inbounds function projector_fourier(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    i_stop = find_truncation_index(psp.β[l][n], tol)
    return hankel_transform(OrbitalLike(), l, psp.r, psp.dr, psp.β[l][n]; i_stop)
end

@inbounds function pseudo_orbital_fourier(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    isnothing(psp.ϕ̃) && return _ -> nothing
    i_stop = find_truncation_index(psp.ϕ̃[l][n], tol)
    return hankel_transform(OrbitalLike(), l, psp.r, psp.dr, psp.ϕ̃[l][n]; i_stop)
end

function valence_charge_density_fourier(psp::NumericPsP; tol=nothing)
    i_stop = find_truncation_index(psp.ρval, tol)
    return hankel_transform(DensityLike(), psp.r, psp.dr, psp.ρval; i_stop)
end

function core_charge_density_fourier(psp::NumericPsP; tol=nothing)
    i_stop = find_truncation_index(psp.ρcore, tol)
    return hankel_transform(DensityLike(), psp.r, psp.dr, psp.ρcore; i_stop)
end

@inbounds function pseudo_energy_correction(T::Type, psp::NumericPsP; tol=nothing)
    i_start = firstindex(psp.Vloc)
    i_stop = find_truncation_index(psp.Vloc, tol)
    integrand(i::Int) = psp.r[i] * (psp.r[i] * psp.Vloc[i] + psp.Zval)
    return T(4π * dotprod(integrand, i_start, i_stop, psp.dr))
end
