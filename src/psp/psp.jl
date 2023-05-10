@doc raw"""
Abstract type representing a pseudopotential.

The structure of the data should facilitate efficient computations.

Required methods:

```julia
# A unique string, usually a hash or checksum.
function identifier(file::AbstractPsP)::AbstractString end
# The symbol of the element for which the pseudopotential is constructed (e.g. `"Ag"`)
function elemental_symbol(psp::AbstractPsP)::AbstractString end
# The maximum angular momentum channel
function max_angular_momentum(psp::AbstractPsP)::Integer end
# The number of non-local projector radial parts for angular momentum `l`
function n_beta_projector_radials(psp::AbstractPsP, l::Integer)::Integer end
# The number of chi function radial parts for angular momentum `l`
function n_chi_projector_radials(psp::AbstractPsP, l::Integer)::Integer end
# The pseudo-atomic valence charge
function valence_charge(psp::AbstractPsP)::Real end
# The charge of the atom which was pseudized (e.g. 8 for Oxygen)
function atomic_charge(psp::AbstractPsP)::Real end
# Whether the pseudopotential is a norm-conserving pseudopotential
function is_norm_conserving(psp::AbstractPsP)::Bool end
# Whether the pseudopotential is an ultrasoft pseudopotential
function is_ultrasoft(psp::AbstractPsP)::Bool end
# Whether the pseudopotential is a projector-augmented wave pseudopotential
function is_paw(psp::AbstractPsP)::Bool end
# Whether the pseudopotential supports spin-orbit coupled calculations
function has_spin_orbit(psp::AbstractPsP)::Bool end
# Whether the pseudopotential contains a core charge density (i.e. supports non-linear core
# correction)
function has_core_density(psp::AbstractPsP)::Bool end
# Whether the pseudopotential contains a valence charge density (i.e. has support for
# constructing a tailored guess charge density)
function has_valence_density(psp::AbstractPsP)::Bool end
# Whether pseudopotential contains chi functions for the valence electrons (i.e.
# has support for computing tailored orbital-projected quantitites)
function has_chi_projectors(psp:AbstractPsP)::Bool end
# The projector coupling coefficients for angular momentum `l`
function beta_projector_coupling(psp::AbstractPsP, l::Integer)::Matrix{Real} end
# Radial distance where the local potential decays to zero within a tolerance `tol`
function local_potential_cutoff_radius(psp::AbstractPsP; tol) end
# Radial distance where the `n`th non-local projector at angular momentum `l` decays to
# zero within a tolerance `tol`
function beta_projector_cutoff_radius(psp::AbstractPsP, l, n; tol) end
# Radial distance where the `n`th chi function at angular momentum `l` decays to
# zero within a tolerance `tol`
function chi_projector_cutoff_radius(psp::AbstractPsP, l, n; tol) end
# Radial distance where the valence charge density decays to zero within a tolerance `tol`
function valence_charge_density_cutoff_radius(psp::AbstractPsP; tol) end
# Radial distance where the core charge density decays to zero within a tolerance `tol`
function core_charge_density_cutoff_radius(psp::AbstractPsP; tol) end
# Returns a function which evaulates the local potential at a real-space radial coordinate
function local_potential_real(psp::AbstractPsP) end
# Returns a function which evaulates the `n`th non-local projector with angular momentum
# `l` at a real-space radial coordinate
function beta_projector_real(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the `n`th chi function with angular momentum
# `l` at a real-space radial coordinate
function chi_projector_real(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the valence charge density at a real-space
# radial coordinate
function valence_charge_density_real(psp::AbstractPsP) end
# Returns a function which evaulates the core charge density at a real-space
# radial coordinate
function core_charge_density_real(psp::AbstractPsP) end
# Returns a function which evaulates the local potential at a fourier-space radial
# coordinate
function local_potential_fourier(psp::AbstractPsP) end
# Returns a function which evaulates the `n`th non-local projector with angular momentum
# `l` at a fourier-space radial coordinate
function beta_projector_fourier(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the `n`th chi function with angular momentum
# `l` at a fourier-space radial coordinate
function chi_projector_fourier(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the valence charge density at a fourier-space
# radial coordinate
function valence_charge_density_fourier(psp::AbstractPsP) end
# Returns a function which evaulates the core charge density at a fourier-space
# radial coordinate
function core_charge_density_fourier(psp::AbstractPsP) end
# The pseudo-potential energy correction
function psp_energy_correction(psp::AbstractPsP) end
```
"""
abstract type AbstractPsP end

#!!! Required functions !!!#
"""
Identifying data (preferably unique).
"""
function identifier(psp::AbstractPsP) end

"""
Element which the pseudopotential was constructed to reproduce.
"""
function element(psp::AbstractPsP) end

"""
Maximum angular momentum channel of the pseudopotential.
"""
function max_angular_momentum(psp::AbstractPsP) end

"""
Number of radial functions of a given quantity provided by the pseudopotential at angular
momentum ``l``.
"""
function n_radials(psp::AbstractPsP, quantity::AbstractPsPQuantity, l) end

"""
Pseudo-atomic valence charge.
"""
function valence_charge(psp::AbstractPsP) end

"""
Charge of the atom corresponding to the `psp`'s `element`.
"""
function atomic_charge(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(psp::AbstractPsP) end

"""
Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(psp::AbstractPsP) end

"""
Find the cutoff radius of a given quantity.
The optional absolute tolerance `tol` can be used to find the radial coordinate where the
absolute value of the quantity decays to below the tolerance rather than returning the last
index of a vector or a parsed cutoff radius.

If no quantity is provided, the keyword argument `f` can be used to select whether the
minimum or maximum cutoff radius over all pseudopotential quantities is returned.
"""
function cutoff_radius(psp::AbstractPsP, quantity::AbstractPsPQuantity; tol=nothing) end
function cutoff_radius(psp::AbstractPsP, quantity::PsPProjector, l, n; tol=nothing) end

"""
Construct a callable which evaluates the requested pseudopotential quantity in either
Fourier- or real-space.
"""
function psp_quantity_evaluator(psp::AbstractPsP, quantity::AbstractPsPQuantity,
                                space::EvaluationSpace) end
function psp_quantity_evaluator(psp::AbstractPsP, quantity::AbstractPsPQuantity, l, n,
                                space::EvaluationSpace) end

"""
Pseudopotential energy correction (the DC component of the Fourier transform of the
local part of the pseudopotential).
"""
function psp_energy_correction(T::Type, psp::AbstractPsP) end

@doc raw"""
Get the internal representation of a given quantity, if available.
"""
function get_quantity(psp::AbstractPsP, quantity::AbstractPsPQuantity) end
function get_quantity(psp::AbstractPsP, quantity::PsPProjector, l, n) end

"""
Check that a pseudopotential has a given quantity.
"""
function has_quantity(psp::AbstractPsP, quantity::AbstractPsPQuantity) end

#!!! Convenience functions !!!#
"""
Angular momenta values of the pseudopotential.
"""
angular_momenta(psp::AbstractPsP) = 0:max_angular_momentum(psp)

"""
Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(psp::AbstractPsP)::Symbol = has_spin_orbit(psp) ? :full : :scalar

"""
Formalism of the pseudopotential (norm-conserving, ultrasoft, projector-augmented wave,
or Coulomb).
"""
function formalism(psp::AbstractPsP)::Type
    is_norm_conserving(psp) && return NormConservingPsP
    is_ultrasoft(psp) && return UltrasoftPsP
    is_paw(psp) && return ProjectorAugmentedWavePsP
end

function n_radials(psp::AbstractPsP, quantity::PsPProjector)
    return sum(l -> n_radials(psp, quantity, l), angular_momenta(psp); init=0)
end

"""
The number of radial + angular parts corresponding to a quantity at a given angular
momentum, or the total number of radial + angular parts for all angular momenta for
a given quantity.

i.e., count the number of combinations for the quantum numbers ``l`` and ``m`` up to the
maximum angular momentum of the pseudopotential, accounting for multi-projector
pseudopotentials that provide multiple radial parts at an individual angular momentum.
"""
function n_angulars(psp::AbstractPsP, quantity::PsPProjector, l)
    return (2l + 1) * n_radials(psp, quantity, l)
end
function n_angulars(psp::AbstractPsP, quantity::PsPProjector)
    return sum(l -> n_angulars(psp, quantity, l), angular_momenta(psp); init=0)
end

function cutoff_radius(psp::AbstractPsP, quantity::PsPProjector, l::Int; f=minimum,
                       tol=nothing)
    cutoff_radii = map(n -> cutoff_radius(psp, quantity, l, n; tol),
                       1:n_radials(psp, quantity, l))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    return isempty(cutoff_radii) ? nothing : f(cutoff_radii)
end

function cutoff_radius(psp::AbstractPsP, quantity::PsPProjector; f=minimum, tol=nothing)
    cutoff_radii = map(l -> cutoff_radius(psp, quantity, l; tol), angular_momenta(psp))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    return isempty(cutoff_radii) ? nothing : f(cutoff_radii)
end

function cutoff_radius(psp::AbstractPsP,
                       quantities::AbstractVector{AbstractPsPQuantity}=[ValenceChargeDensity(),
                                                                        CoreChargeDensity(),
                                                                        BetaProjector(),
                                                                        ChiProjector(),
                                                                        LocalPotential()];
                       f=minimum, tol=nothing)
    cutoff_radii = map(q -> cutoff_radius(psp, q; tol), quantities)
    cutoff_radii = filter(!isnothing, cutoff_radii)
    return isempty(cutoff_radii) ? nothing : f(cutoff_radii)
end

Base.Broadcast.broadcastable(psp::AbstractPsP) = Ref(psp)

# TODO: figure out a better equality metric
# Base.:(==)(psp1::AbstractPsP, psp2::AbstractPsP) = identifier(psp1) == identifier(psp2)
# Base.hash(psp::AbstractPsP) = hash(psp.checksum)

function Base.show(io::IO, psp::AbstractPsP)
    typename = string(typeof(psp))
    el = element(psp)
    z = valence_charge(psp)
    return print(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", psp::AbstractPsP)
    println(io, typeof(psp))
    @printf "%032s: %s\n" "identifier" identifier(psp)
    @printf "%032s: %s\n" "formalism" formalism(psp)
    @printf "%032s: %s\n" "element" element(psp)
    @printf "%032s: %f\n" "valence charge" valence_charge(psp)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(psp)
    @printf "%032s: %s\n" "non-linear core correction" has_quantity(psp,
                                                                    CoreChargeDensity())
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(psp)
    @printf "%032s: %s\n" "number of beta projectors" n_radials(psp, BetaProjector())
    @printf "%032s: %s" "number of chi projectors" n_radials(psp, ChiProjector())
end
