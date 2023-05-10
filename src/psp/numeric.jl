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
function n_radials(psp::NumericPsP, quantity::PsPProjector, l)
    return has_quantity(psp, quantity) ? length(get_quantity(psp, quantity, l)) : 0
end

get_quantity(psp::NumericPsP, ::LocalPotential) = psp.Vloc
function get_quantity(psp::NumericPsP, quantity::PsPProjector, l, n)
    return get_quantity(psp, quantity, l)[n]
end
get_quantity(psp::NumericPsP, quantity::PsPProjector, l) = get_quantity(psp, quantity)[l]
get_quantity(psp::NumericPsP, ::BetaProjector) = psp.β
get_quantity(psp::NumericPsP, ::BetaCoupling) = psp.D
get_quantity(psp::NumericPsP, ::BetaCoupling, l) = psp.D[l]
get_quantity(psp::NumericPsP, ::BetaCoupling, l, n) = psp.D[l][n, n]
get_quantity(psp::NumericPsP, ::BetaCoupling, l, n, m) = psp.D[l][n, m]
get_quantity(psp::NumericPsP, ::ChiProjector) = psp.χ
get_quantity(psp::NumericPsP, ::ValenceChargeDensity) = psp.ρval
get_quantity(psp::NumericPsP, ::CoreChargeDensity) = psp.ρcore
get_quantity(::NumericPsP, ::AugmentationFunction) = nothing

function has_quantity(psp::NumericPsP, quantity::AbstractPsPQuantity)
    if !isnothing(get_quantity(psp, quantity))  # First check if the quantity is explicit nothing
        !isempty(get_quantity(psp, quantity)) && return true  # Then check if it is empty
    end
    return false
end

function cutoff_radius(psp::NumericPsP, quantity::AbstractPsPQuantity; f=nothing,
                       tol=nothing)
    !has_quantity(psp, quantity) && return nothing
    f = get_quantity(psp, quantity)
    return psp.r[find_truncation_index(f, tol)]
end

function cutoff_radius(psp::NumericPsP, quantity::PsPProjector, l, n; f=nothing,
                       tol=nothing)
    !has_quantity(psp, quantity) && return nothing
    f = get_quantity(psp, quantity, l, n)
    return psp.r[find_truncation_index(f, tol)]
end

function cutoff_radius(psp::NumericPsP, quantity::PsPProjector; f=minimum, tol=nothing)
    cutoff_radii = map(l -> cutoff_radius(psp, quantity, l; tol), angular_momenta(psp))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    return isempty(cutoff_radii) ? nothing : f(cutoff_radii)
end

# Real-space angular-independent quantities: local potential, atomic charge densities
function psp_quantity_evaluator(psp::NumericPsP, quantity::AbstractPsPQuantity, ::RealSpace)
    !has_quantity(psp, quantity) && return _ -> nothing
    return build_interpolator_real(get_quantity(psp, quantity), psp.r)
end

# Real-space angular-dependent quantities: projectors
function psp_quantity_evaluator(psp::NumericPsP, quantity::PsPProjector, l, n, ::RealSpace)
    !has_quantity(psp, quantity) && return _ -> nothing
    return build_interpolator_real(get_quantity(psp, quantity, l, n), psp.r)
end

# TODO: r, dr, and Δr need to be truncated to match the length of the underlying
# TODO: quantity as stored in the psp
# Fourier-space angular-independent quantities: atomic charge densities
function psp_quantity_evaluator(psp::NumericPsP, quantity::AbstractPsPQuantity, space::FourierSpace; quadrature_method=Simpson())
    !has_quantity(psp, quantity) && return _ -> nothing
    istop = lastindex(get_quantity(psp, quantity))
    r = @view psp.r[begin:istop]
    dr = @view psp.dr[begin:istop]
    if typeof(psp.Δr) <: AbstractVector
        Δr = @view psp.Δr[begin:istop-1]
    else
        Δr = psp.Δr
    end
    return psp_quantity_evaluator(psp, quantity, space, r, dr, Δr; quadrature_method)
end

function psp_quantity_evaluator(psp::NumericPsP, quantity::AbstractPsPQuantity,
                                ::FourierSpace, r::AbstractVector,
                                dr::AbstractVector,
                                Δr::Union{Real,AbstractVector};
                                quadrature_method=Simpson())
    !has_quantity(psp, quantity) && return _ -> nothing
    f = psp_quantity_evaluator(psp, quantity, RealSpace())
    # l=0 for angular-independent quantities
    return hankel_transform(r, f, dr, Δr, 0; quadrature_method)
end

# Fourier-space angular-dependent quantities: projectors
function psp_quantity_evaluator(psp::NumericPsP, quantity::PsPProjector, l, n, space::FourierSpace; quadrature_method=Simpson())
    !has_quantity(psp, quantity) && return _ -> nothing
    istop = lastindex(get_quantity(psp, quantity, l, n))
    r = @view psp.r[begin:istop]
    dr = @view psp.dr[begin:istop]
    if typeof(psp.Δr) <: AbstractVector
        Δr = @view psp.Δr[begin:istop-1]
    else
        Δr = psp.Δr
    end
    return psp_quantity_evaluator(psp, quantity, l, n, space, r, dr, Δr; quadrature_method)
end

function psp_quantity_evaluator(psp::NumericPsP, quantity::PsPProjector, l, n,
                                ::FourierSpace, r::AbstractVector,
                                dr::AbstractVector,
                                Δr::Union{Real,AbstractVector};
                                quadrature_method=Simpson())
    !has_quantity(psp, quantity) && return _ -> nothing
    f = psp_quantity_evaluator(psp, quantity, l, n, RealSpace())
    return hankel_transform(r, f, dr, Δr, l; quadrature_method)
end

function psp_quantity_evaluator(psp::NumericPsP, quantity::LocalPotential, space::FourierSpace; quadrature_method=Simpson(), correction_method=CoulombCorrection())
    istop = lastindex(get_quantity(psp, quantity))
    r = @view psp.r[begin:istop]
    dr = @view psp.dr[begin:istop]
    if typeof(psp.Δr) <: AbstractVector
        Δr = @view psp.Δr[begin:istop-1]
    else
        Δr = psp.Δr
    end
    return psp_quantity_evaluator(psp, quantity, space, r, dr, Δr; quadrature_method, correction_method)
end

@inbounds function psp_quantity_evaluator(psp::NumericPsP, ::LocalPotential, ::FourierSpace,
                                          r::AbstractVector,
                                          dr::AbstractVector,
                                          Δr::Union{Real,AbstractVector};
                                          quadrature_method=Simpson(),
                                          correction_method=CoulombCorrection())
    function Vloc(q::T) where {T}
        # For q > 0, compute the Hankel-Fourier transform of the local potential
        if !iszero(q)
            # Pre-generate an evaluation function for the local potential in real space
            Vloc = psp_quantity_evaluator(psp, LocalPotential(), RealSpace())
            # Pre-generate the real-space correction function with a single argument
            corr(r) = local_potential_correction(psp, correction_method, RealSpace(), r)
            # Define the integrand for the Hankel transform, subtracting a Coulomb-like
            # correction term which should remove a 1/r-like tail from the local potential
            # and localize it well enough that the integral is well-behaved
            # f(r) = Vloc(r) - c(r) / r
            function integrand(r::Real)
                # Note that Vloc(r) is stored _without_ an r² prefactor, unlike all other
                # radial quantities!
                # Hankel transform   : 4π ∫ r² f(r) jₗ(q * r) dr
                # Correction         : c(r) / r
                # jₗ(q * r) for l=0   : sin(q * r) / (q * r)
                # Substitute         : 4π ∫ r² (Vloc(r) - c(r) / r) sin(q * r) / (q * r) dr
                # Cancel and reorder : 4π / q ∫ (r * Vloc(r) - c(r)) sin(q * r) dr
                # Integrand          : (r * Vloc(r) - c(r)) sin(q * r)
                return (r * Vloc(r) - corr(r)) * sin(q * r)
            end
            integral = integrate(r, integrand, dr, Δr, quadrature_method)
            # Compute the analytical Fourier-Hankel transform of the correction
            corr = 4T(π) *
                   local_potential_correction(psp, correction_method, FourierSpace(), q)
            # Recall that 4π / q is factored out of the integral above
            return 4T(π) / q * integral + corr
        end
        # For q = 0, return the pseudo. energy correction (to avoid division by zero)
        return psp_energy_correction(T, psp; quadrature_method, correction_method)
    end
    Vloc(Q::AbstractVector) = Vloc(norm(Q))
    return Vloc
end

function local_potential_correction(psp::NumericPsP, ::ErfCorrection, ::RealSpace, r::Real)
    return -psp.Zval * erf(r)  # As in QuantumESPRESSO
end
function local_potential_correction(psp::NumericPsP, ::ErfCorrection, ::FourierSpace,
                                    q::Real)
    return -psp.Zval * exp(-q^2 / 4) / q^2  # As in QuantumESPRESSO
end

function local_potential_correction(psp::NumericPsP, ::CoulombCorrection, ::RealSpace,
                                    ::Real)
    return -psp.Zval  # As in ABINIT
end
function local_potential_correction(psp::NumericPsP, ::CoulombCorrection, ::FourierSpace,
                                    q::Real)
    return -psp.Zval / q^2  # As in ABINIT
end

function psp_energy_correction(T::Type{<:Real}, psp::NumericPsP; quadrature_method=Simpson(), correction_method=CoulombCorrection())
    istop = lastindex(get_quantity(psp, LocalPotential()))
    r = @view psp.r[begin:istop]
    dr = @view psp.dr[begin:istop]
    if typeof(psp.Δr) <: AbstractVector
        Δr = @view psp.Δr[begin:istop-1]
    else
        Δr = psp.Δr
    end
    return psp_energy_correction(T, psp, r, dr, Δr; quadrature_method, correction_method)
end

@inbounds function psp_energy_correction(T::Type{<:Real}, psp::NumericPsP,
                                         r::AbstractVector,
                                         dr::AbstractVector,
                                         Δr::Union{Real,AbstractVector};
                                         quadrature_method=Simpson(),
                                         correction_method=CoulombCorrection())
    Vloc = psp_quantity_evaluator(psp, LocalPotential(), RealSpace())
    corr(r) = local_potential_correction(psp, correction_method, RealSpace(), r)
    return 4T(π) *
           integrate(r, ri -> ri * (ri * Vloc(ri) - corr(ri)), dr, Δr, quadrature_method)
end
