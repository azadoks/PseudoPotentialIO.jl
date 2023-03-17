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
# Checksum
checksum::Vector{UInt8}
# Atomic total charge in units of electron charge
Zatom::Number
# Pseudo-atomic valence charge in units of electron charge
Zval::Number
# Maximum angular momentum
lmax::Integer
# Radial mesh in units of Bohr
r::AbstractVector{Real}
# Radial mesh spacing in units of Bohr
dr::Union{Real, AbstractVector{Real}}
# Local potential on the radial mesh in units of Hartree (without r² prefactor)
Vloc::AbstractVector{Real}
## The units of `D` and `β` should be such that `⟨ βˡₙ | Dˡₙₙ | βˡₙ ⟩` gives Hartree
# Nonlocal projector coupling constants D[l][n,n']
D::OffsetVector{AbstractMatrix{Real}}
# Nonlocal projectors on the radial mesh, multiplied by the mesh squared: r²β[l][n]
β::OffsetVector{AbstractVector{AbstractVector{Real}}}

## "Optional" fields (must still exist, but could be Union{Nothing})
# Model core charge density (non-linear core correction) on the radial mesh, multiplied by
# the mesh squared: r²ρcore
ρcore::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic valence charge density on the radial mesh, multiplied by the mesh squared:
# r²ρval
ρval::Union{Nothing,AbstractVector{Real}}
# Pseudo-atomic orbitals on the radial mesh, multiplied by the mesh squared: r²χ[l][n]
χ::Union{Nothing,OffsetVector{AbstractVector{AbstractVector{Real}}}}
```
"""
abstract type NumericPsP{T} <: AbstractPsP end

identifier(psp::NumericPsP)::String = bytes2hex(psp.checksum)
function element(psp::NumericPsP)::String
    return PeriodicTable.elements[Int(psp.Zatom)].symbol
end
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
n_projector_radials(psp::NumericPsP, l::Int)::Int = length(psp.β[l])
n_chi_function_radials(psp::NumericPsP, l::Int)::Int = isnothing(psp.χ) ? 0 :
                                                         length(psp.χ[l])
valence_charge(psp::NumericPsP) = psp.Zval
atomic_charge(psp::NumericPsP) = psp.Zatom
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation
has_core_density(psp::NumericPsP)::Bool = !isnothing(psp.ρcore)
has_valence_density(psp::NumericPsP)::Bool = !isnothing(psp.ρval)
has_chi_functions(psp::NumericPsP)::Bool = !isnothing(psp.χ)

function local_potential_cutoff_radius(psp::NumericPsP; tol=nothing)
    ir_cut = find_truncation_index(psp.Vloc, tol)
    return psp.r[ir_cut]
end

function projector_cutoff_radius(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    ir_cut = find_truncation_index(psp.β[l][n], tol)
    return psp.r[ir_cut]
end

function chi_function_cutoff_radius(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    !has_chi_functions(psp) && return nothing
    ir_cut = find_truncation_index(psp.χ[l][n], tol)
    return psp.r[ir_cut]
end

function valence_charge_density_cutoff_radius(psp::NumericPsP; tol=nothing)
    !has_valence_density(psp) && return nothing
    ir_cut = find_truncation_index(psp.ρval, tol)
    return psp.r[ir_cut]
end

function core_charge_density_cutoff_radius(psp::NumericPsP; tol=nothing)
    !has_core_density(psp) && return nothing
    ir_cut = find_truncation_index(psp.ρcore, tol)
    return psp.r[ir_cut]
end

function projector_coupling(psp::NumericPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.D[l]
end

function local_potential_real(psp::NumericPsP)
    return build_interpolator_real(psp.Vloc, psp.r)
end

function projector_real(psp::NumericPsP, l::Int, n::Int)
    return build_interpolator_real(psp.β[l][n], psp.r)
end

function chi_function_real(psp::NumericPsP, l::Int, n::Int)
    isnothing(psp.χ) && return _ -> nothing
    return build_interpolator_real(psp.χ[l][n], psp.r)
end

function valence_charge_density_real(psp::NumericPsP)
    isnothing(psp.ρval) && return _ -> nothing
    return build_interpolator_real(psp.ρval, psp.r)
end

function core_charge_density_real(psp::NumericPsP)
    isnothing(psp.ρcore) && return _ -> nothing
    return build_interpolator_real(psp.ρcore, psp.r)
end

@inbounds function local_potential_fourier(psp::NumericPsP; tol=nothing)
    i_start = firstindex(psp.Vloc)
    i_stop = find_truncation_index(psp.Vloc, tol)
    function Vloc(q)
        function integrand(i::Int)
            return psp.r[i] * fast_sphericalbesselj0(q * psp.r[i]) *
                   (psp.r[i] * psp.Vloc[i] + psp.Zval)
        end
        integral = dotprod(integrand, i_start, i_stop, psp.dr)
        return 4π * (integral - psp.Zval / q^2)
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

@inbounds function projector_fourier(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    i_stop = find_truncation_index(psp.β[l][n], tol)
    return hankel_transform(psp.β[l][n], l, psp.r, psp.dr; i_stop)
end

@inbounds function chi_function_fourier(psp::NumericPsP, l::Int, n::Int; tol=nothing)
    isnothing(psp.χ) && return _ -> nothing
    i_stop = find_truncation_index(psp.χ[l][n], tol)
    return hankel_transform(psp.χ[l][n], l, psp.r, psp.dr; i_stop)
end

function valence_charge_density_fourier(psp::NumericPsP; tol=nothing)
    i_stop = find_truncation_index(psp.ρval, tol)
    return hankel_transform(psp.ρval, 0, psp.r, psp.dr; i_stop)
end

function core_charge_density_fourier(psp::NumericPsP; tol=nothing)
    i_stop = find_truncation_index(psp.ρcore, tol)
    return hankel_transform(psp.ρcore, 0, psp.r, psp.dr; i_stop)
end

@inbounds function pseudo_energy_correction(T::Type, psp::NumericPsP; tol=nothing)
    i_start = firstindex(psp.Vloc)
    i_stop = find_truncation_index(psp.Vloc, tol)
    integrand(i::Int) = psp.r[i] * (psp.r[i] * psp.Vloc[i] + psp.Zval)
    return T(4π * dotprod(integrand, i_start, i_stop, psp.dr))
end
