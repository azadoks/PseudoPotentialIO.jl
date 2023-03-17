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
function n_projector_radials(psp::AbstractPsP, l::Integer)::Integer end
# The number of chi function radial parts for angular momentum `l`
function n_chi_function_radials(psp::AbstractPsP, l::Integer)::Integer end
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
function has_chi_functions(psp:AbstractPsP)::Bool end
# The projector coupling coefficients for angular momentum `l`
function projector_coupling(psp::AbstractPsP, l::Integer)::Matrix{Real} end
# Radial distance where the local potential decays to zero within a tolerance `tol`
function local_potential_cutoff_radius(psp::AbstractPsP; tol) end
# Radial distance where the `n`th non-local projector at angular momentum `l` decays to
# zero within a tolerance `tol`
function projector_cutoff_radius(psp::AbstractPsP, l, n; tol) end
# Radial distance where the `n`th chi function at angular momentum `l` decays to
# zero within a tolerance `tol`
function chi_function_cutoff_radius(psp::AbstractPsP, l, n; tol) end
# Radial distance where the valence charge density decays to zero within a tolerance `tol`
function valence_charge_density_cutoff_radius(psp::AbstractPsP; tol) end
# Radial distance where the core charge density decays to zero within a tolerance `tol`
function core_charge_density_cutoff_radius(psp::AbstractPsP; tol) end
# Returns a function which evaulates the local potential at a real-space radial coordinate
function local_potential_real(psp::AbstractPsP) end
# Returns a function which evaulates the `n`th non-local projector with angular momentum
# `l` at a real-space radial coordinate
function projector_real(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the `n`th chi function with angular momentum
# `l` at a real-space radial coordinate
function chi_function_real(psp::AbstractPsP, l, n) end
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
function projector_fourier(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the `n`th chi function with angular momentum
# `l` at a fourier-space radial coordinate
function chi_function_fourier(psp::AbstractPsP, l, n) end
# Returns a function which evaulates the valence charge density at a fourier-space
# radial coordinate
function valence_charge_density_fourier(psp::AbstractPsP) end
# Returns a function which evaulates the core charge density at a fourier-space
# radial coordinate
function core_charge_density_fourier(psp::AbstractPsP) end
# The pseudo-potential energy correction
function pseudo_energy_correction(psp::AbstractPsP) end
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
Number of radial parts `Rₗₙ(|r|)` of the Kleinman-Bylander projectors `Rₗₙ(|r|)Yₗₘ(r̂)` at a
given angular momentum `l`.
"""
function n_projector_radials(psp::AbstractPsP, l) end

"""
Number of radial parts of the chi-functions with angular momentum `l`.
"""
function n_chi_function_radials(psp::AbstractPsP, l) end

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
Whether the pseudopotential contains non-linear core correction data (model core charge
density).
"""
function has_core_density(psp::AbstractPsP) end

"""
Whether the pseudopotential contains valence charge density data.
"""
function has_valence_density(psp::AbstractPsP) end

"""
Whether the pseudopotential contains pseudoatomic orbitals.
"""
function has_chi_functions(psp::AbstractPsP) end

"""
Projector coupling matrix at angular momentum `l`.
"""
function projector_coupling(psp::AbstractPsP, l) end

"""
Cutoff radius of the local potential in real-space.
"""
function local_potential_cutoff_radius(psp::AbstractPsP; tol=nothing) end

"""
Cutoff radius of the `n`th Kleinman-Bylander non-local projector at angular momentum `l` in
real-space.
"""
function projector_cutoff_radius(psp::AbstractPsP, l::Int, n::Int; tol=nothing) end

"""
Cutoff radius of the `n`th chi function at angular momentum `l` in real-space.
"""
function chi_function_cutoff_radius(psp::AbstractPsP, l::Int, n::Int; tol=nothing) end

"""
Cutoff radius of the valence charge density in real-space.
"""
function valence_charge_density_cutoff_radius(psp::AbstractPsP; tol=nothing) end

"""
Cutoff radius of the core charge density in real-space.
"""
function core_charge_density_cutoff_radius(psp::AbstractPsP; tol=nothing) end

"""
Local part of the pseudopotential evaluated at real-space point `r`.
"""
function local_potential_real(psp::AbstractPsP)
    return _ -> nothing
end

"""
The `n`th nonlocal Kleinman-Bylander projector at angular momentum `l` evaluated at
real-space point `r`.
"""
function projector_real(psp::AbstractPsP, l::Integer, n::Integer)
    return _ -> nothing
end

"""
The `n`th chi function at angular momentum `l` evaulated at real-space point `r`.
"""
function chi_function_real(psp::AbstractPsP, l::Integer, n::Integer)
    return _ -> nothing
end

"""
Pseudo-atomic valence charge density evaluated at real-space point `r`.
"""
function valence_charge_density_real(psp::AbstractPsP)
    return _ -> nothing
end

"""
Model core charge density evaluated at real-space point `r`.
"""
function core_charge_density_real(psp::AbstractPsP)
    return _ -> nothing
end

"""
Local part of the pseudopotential evaluated at reciprocal-space point `q`.
"""
function local_potential_fourier(psp::AbstractPsP)
    return _ -> nothing
end

"""
The `n`th nonlocal Kleinman-Bylander projector at angular momentum `l` evaluated at
reciprocal-space point `q`.
"""
function projector_fourier(psp::AbstractPsP, l::Integer, n::Integer)
    return _ -> nothing
end

"""
The `n`th chi function at angular momentum `l` evaulated at reciprocal-space
point `q`.
"""
function chi_function_fourier(psp::AbstractPsP, l::Integer, n::Integer)
    return _ -> nothing
end

"""
Pseudo-atomic valence charge density evaluated at reciprocal-space point `q`.
"""
function valence_charge_density_fourier(psp::AbstractPsP)
    return _ -> nothing
end

"""
Model core charge density evaluated at reciprocal-space point `q`.
"""
function core_charge_density_fourier(psp::AbstractPsP)
    return _ -> nothing
end

"""
Pseudopotential energy correction (the DC component of the Fourier transform of the
local part of the pseudopotential).
"""
function pseudo_energy_correction(T::Type, psp::AbstractPsP) end

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

"""
Indices of projector radial parts at a given angular momentum.
"""
projector_radial_indices(psp::AbstractPsP, l) = 1:n_projector_radials(psp, l)

"""
Number of radial parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_projector_radials(psp::AbstractPsP)
    return sum(l -> n_projector_radials(psp, l), angular_momenta(psp); init=0)
end

"""
Number of angular parts `Yₗₘ(r̂)` of the Kleinman-Bylander projectors `Rₗₙ(|r|)Yₗₘ(r̂)` at a
given angular momentum `l`.
"""
function n_projector_angulars(psp::AbstractPsP, l)
    # for angular momentum l, magnetic q.n. m ∈ -l:l
    return n_projector_radials(psp::AbstractPsP, l) * (2l + 1)
end

"""
Number of angular parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_projector_angulars(psp::AbstractPsP)
    return sum(l -> n_projector_angulars(psp, l), angular_momenta(psp); init=0)
end

"""
Indices of pseudo-atomic wavefunction radial parts at a given angular momentum.
"""
chi_function_radial_indices(psp::AbstractPsP, l) = 1:n_chi_function_radials(psp, l)

"""
Number chi-functions `Rₗₙ(|r|) * Yₗₘ(r̂)` at angular momenta `l` up to the
maximum angular momentum channel.
"""
function n_chi_function_radials(psp::AbstractPsP)
    return sum(l -> n_chi_function_radials(psp, l), angular_momenta(psp); init=0)
end

"""
Number of angular parts of the chi-functions with angular momentum `l`.
"""
function n_chi_function_angulars(psp::AbstractPsP, l)
    # for angular momentum l, magnetic q.n. m ∈ -l:l
    return n_chi_function_radials(psp::AbstractPsP, l) * (2l + 1)
end

"""
Number of angular parts of the chi-functions at all angular momenta up
to the maximum angular momentum channel.
"""
function n_chi_function_angulars(psp::AbstractPsP)
    return sum(l -> n_pseudo_orbtial_angulars(psp, l), angular_momenta(psp); init=0)
end

function projector_cutoff_radius(psp::AbstractPsP, l::Int; f=minimum, tol=nothing)
    cutoff_radii = map(n -> projector_cutoff_radius(psp, l, n; tol),
                       projector_radial_indices(psp, l))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    isempty(cutoff_radii) && return nothing
    return f(cutoff_radii)
end

function projector_cutoff_radius(psp::AbstractPsP; f=minimum, tol=nothing)
    cutoff_radii = map(l -> projector_cutoff_radius(psp, l; tol),
                       angular_momenta(psp))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    isempty(cutoff_radii) && return nothing
    return f(cutoff_radii)
end

function chi_function_cutoff_radius(psp::AbstractPsP, l::Int; f=minimum, tol=nothing)
    cutoff_radii = map(n -> chi_function_cutoff_radius(psp, l, n; tol),
                       chi_function_radial_indices(psp, l))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    isempty(cutoff_radii) && return nothing
    return f(cutoff_radii)
end

function chi_function_cutoff_radius(psp::AbstractPsP; f=minimum, tol=nothing)
    cutoff_radii = map(l -> chi_function_cutoff_radius(psp, l; tol),
                       angular_momenta(psp))
    cutoff_radii = filter(!isnothing, cutoff_radii)
    isempty(cutoff_radii) && return nothing
    return f(cutoff_radii)
end

"""
Find the cutoff radius for the pseudopotential. Supply a function `f` to determine what kind
of reduction over cutoff radii for different quantities is performed. Supply `tol` if
you'd like to truncate the quantities where they decay to `|f| < tol`.
"""
function pseudo_cutoff_radius(psp::AbstractPsP; f=minimum, tol=nothing)
    cutoff_radii = [local_potential_cutoff_radius(psp; tol),
                    projector_cutoff_radius(psp; f, tol),
                    chi_function_cutoff_radius(psp; f, tol),
                    valence_charge_density_cutoff_radius(psp; tol),
                    core_charge_density_cutoff_radius(psp; tol)]
    cutoff_radii = filter(!isnothing, cutoff_radii)
    return f(cutoff_radii)
end

"""
Projector coupling constant between the `n`th and `m`th projector with angular momentum `l`.
"""
function projector_coupling(psp::AbstractPsP, l::Int, n::Int, m::Int)
    return projector_coupling(psp, l)[n, m]
end

"""
Projector coupling constant between the `n`th projector with angular momentum `l` itself.
"""
function projector_coupling(psp::AbstractPsP, l::Int, n::Int)
    return projector_coupling(psp, l, n, n)
end

Base.Broadcast.broadcastable(psp::AbstractPsP) = Ref(psp)

Base.:(==)(psp1::AbstractPsP, psp2::AbstractPsP) = identifier(psp1) == identifier(psp2)
Base.hash(psp::AbstractPsP) = hash(psp.checksum)

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
    @printf "%032s: %s\n" "non-linear core correction" has_core_density(psp)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(psp)
    @printf "%032s: %s\n" "number of projectors" n_projector_radials(psp)
    @printf "%032s: %s" "number of chi functions" n_chi_function_radials(psp)
end
