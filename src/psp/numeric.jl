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

identifier(psp::NumericPsP)::String = psp.identifier
element(psp::NumericPsP)::PeriodicTable.Element = PeriodicTable.elements[psp.Zatom]
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
valence_charge(psp::NumericPsP) = psp.Zval
atomic_charge(psp::NumericPsP) = psp.Zatom
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation
function n_radials(q::PsPProjector, psp::NumericPsP, l)
    return has_quantity(q, psp) ? length(get_quantity(q, psp, l)) : 0
end

get_quantity(::LocalPotential, psp::NumericPsP) = psp.Vloc
get_quantity(q::PsPProjector, psp::NumericPsP, l, n) = get_quantity(q, psp, l)[n]
get_quantity(q::PsPProjector, psp::NumericPsP, l) = get_quantity(q, psp)[l]
get_quantity(::BetaProjector, psp::NumericPsP) = psp.β
get_quantity(::BetaCoupling, psp::NumericPsP) = psp.D
get_quantity(::BetaCoupling, psp::NumericPsP, l) = psp.D[l]
get_quantity(::BetaCoupling, psp::NumericPsP, l, n) = psp.D[l][n, n]
get_quantity(::BetaCoupling, psp::NumericPsP, l, n, m) = psp.D[l][n, m]
get_quantity(::ChiProjector, psp::NumericPsP) = psp.χ
get_quantity(::ValenceChargeDensity, psp::NumericPsP) = psp.ρval
get_quantity(::CoreChargeDensity, psp::NumericPsP) = psp.ρcore
get_quantity(::AugmentationFunction, psp::NumericPsP) = nothing

has_quantity(q::AbstractPsPQuantity, psp::NumericPsP) = !isnothing(get_quantity(q, psp))

function cutoff_radius(quantity::AbstractPsPQuantity, psp::NumericPsP; f=nothing,
                       tol=nothing)
    !has_quantity(quantity, psp) && return nothing
    f = get_quantity(quantity, psp)
    return psp.r[find_truncation_index(f, tol)]
end

function cutoff_radius(quantity::PsPProjector, psp::NumericPsP, l, n; f=nothing,
                       tol=nothing)
    !has_quantity(quantity, psp) && return nothing
    f = get_quantity(quantity, psp, l, n)
    return psp.r[find_truncation_index(f, tol)]
end

# Real-space angular-independent quantities: local potential, atomic charge densities
function psp_quantity_evaluator(::RealSpace, quantity::AbstractPsPQuantity, psp::NumericPsP)
    !has_quantity(quantity, psp) && return _ -> nothing
    return build_interpolator_real(get_quantity(quantity, psp), psp.r)
end

# Real-space angular-dependent quantities: projectors
function psp_quantity_evaluator(::RealSpace, quantity::PsPProjector, psp::NumericPsP, l, n)
    !has_quantity(quantity, psp) && return _ -> nothing
    return build_interpolator_real(get_quantity(quantity, psp, l, n), psp.r)
end

# Fourier-space angular-independent quantities: atomic charge densities
# Integrate on the internal pseudopotential grid
function psp_quantity_evaluator(::FourierSpace, quantity::AbstractPsPQuantity,
                                psp::NumericPsP; i_start=nothing, i_stop=nothing,
                                integrator=trapezoid)
    !has_quantity(quantity, psp) && return _ -> nothing
    f = get_quantity(quantity, psp)
    i_start = isnothing(i_start) ? firstindex(f) : i_start
    i_stop = isnothing(i_stop) ? lastindex(f) : i_stop
    return hankel_transform(f, 0, psp.r, psp.dr; i_start, i_stop, integrator)
end

# Fourier-space angular-independent quantities: atomic charge densities
# Integrate on an arbitrary grid by evaluating a cubic spline fit of the quantity
function psp_quantity_evaluator(::FourierSpace, quantity::AbstractPsPQuantity,
                                psp::NumericPsP, r, dr; integrator=trapezoid)
    !has_quantity(quantity, psp) && return _ -> nothing
    f_real = psp_quantity_evaluator(RealSpace(), quantity, psp)
    f = f_real.(r)
    return hankel_transform(f, 0, r, dr; integrator)
end

# Fourier-space angular-independent quantities: projectors
# Integrate on the internal pseudopotential grid
function psp_quantity_evaluator(::FourierSpace, quantity::PsPProjector, psp::NumericPsP, l,
                                n; i_start=nothing, i_stop=nothing, integrator=trapezoid)
    !has_quantity(quantity, psp) && return _ -> nothing
    f = get_quantity(quantity, psp, l, n)
    i_start = isnothing(i_start) ? firstindex(f) : i_start
    i_stop = isnothing(i_stop) ? lastindex(f) : i_stop
    return hankel_transform(f, l, psp.r, psp.dr; i_start, i_stop, integrator)
end

# Fourier-space angular-independent quantities: projectors
# Integrate on an arbitrary grid by evaluating a cubic spline fit of the quantity
function psp_quantity_evaluator(::FourierSpace, quantity::PsPProjector, psp::NumericPsP, l,
                                n, r, dr; integrator=trapezoid)
    !has_quantity(quantity, psp) && return _ -> nothing
    f_real = psp_quantity_evaluator(RealSpace(), quantity, psp, l, n)
    f = f_real.(r)
    return hankel_transform(f, l, r, dr; integrator)
end

# TODO: implement q-e style correction (this is ABINIT-style)
@inbounds function psp_quantity_evaluator(::FourierSpace, ::LocalPotential, psp::NumericPsP;
                                          i_start=nothing, i_stop=nothing,
                                          integrator=trapezoid)
    i_start = isnothing(i_start) ? firstindex(psp.Vloc) : i_start
    i_stop = isnothing(i_stop) ? lastindex(psp.Vloc) : i_stop
    function Vloc(q)
        function integrand(i::Int)
            return psp.r[i] * fast_sphericalbesselj0(q * psp.r[i]) *
                   (psp.r[i] * psp.Vloc[i] + psp.Zval)
        end
        integral = integrator(integrand, i_start, i_stop, psp.r, psp.dr)
        return 4π * (integral - psp.Zval / q^2)
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

@inbounds function psp_quantity_evaluator(::FourierSpace, ::LocalPotential, psp::NumericPsP,
                                          r, dr; integrator=trapezoid)
    f_real = psp_quantity_evaluator(RealSpace(), LocalPotential(), psp)
    f = f_real.(r)
    i_start = firstindex(f)
    i_stop = lastindex(f)
    function Vloc(q)
        function integrand(i::Int)
            return r[i] * fast_sphericalbesselj0(q * r[i]) * (r[i] * f[i] + psp.Zval)
        end
        integral = integrator(integrand, i_start, i_stop, r, dr)
        return 4π * (integral - psp.Zval / q^2)
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

@inbounds function psp_energy_correction(T::Type{<:Real}, psp::NumericPsP; i_start=nothing,
                                         i_stop=nothing, integrator=trapezoid)
    i_start = isnothing(i_start) ? firstindex(psp.Vloc) : i_start
    i_stop = isnothing(i_stop) ? lastindex(psp.Vloc) : i_stop
    integrand(i::Int) = psp.r[i] * (psp.r[i] * psp.Vloc[i] + psp.Zval)
    return 4T(π) * integrator(integrand, i_start, i_stop, psp.r, psp.dr)
end

@inbounds function psp_energy_correction(T::Type{<:Real}, psp::NumericPsP, r, dr;
                                         integrator=trapezoid)
    f_real = psp_quantity_evaluator(RealSpace(), LocalPotential(), psp)
    f = f_real.(r)
    i_start = firstindex(f)
    i_stop = lastindex(f)
    integrand(i::Int) = r[i] * (r[i] * f[i] + psp.Zval)
    return 4T(π) * integrator(integrand, i_start, i_stop, r, dr)
end
