@doc raw"""
Abstract type representing a pseudopotential.

The structure of the data should facilitate efficient computations.

Required methods:

```julia
# The symbol of the element for which the pseudopotential is constructed (e.g. `"Ag"`)
function elemental_symbol(psp::AbstractPsP)::AbstractString end
# The maximum angular momentum channel
function max_angular_momentum(psp::AbstractPsP)::Integer end
# The number of non-local projector radial parts for angular momentum `l`
function n_projector_radials(psp::AbstractPsP, l::Integer)::Integer end
# The number of pseudo-atomic orbital radial parts for angular momentum `l`
function n_pseudo_orbital_radials(psp::AbstractPsP, l::Integer)::Integer end
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
# Whether pseudopotential contains pseudo-atomic orbitals for the valence electrons (i.e.
# has support for computing tailored orbital-projected quantitites)
function has_pseudo_orbitals(psp:AbstractPsP)::Bool end
# The projector coupling coefficients for angular momentum `l`
function projector_coupling(psp::AbstractPsP, l::Integer)::Matrix{Real} end
```
"""
abstract type AbstractPsP end

#!!! Required functions !!!#
"""
Element which the pseudopotential was constructed to reproduce.
"""
function element(psp::AbstractPsP) end

"""
Maximum angular momentum channel in the local part of the pseudopotential.
"""
function max_angular_momentum(psp::AbstractPsP) end

"""
Number of radial parts `Rₗₙ(|r|)` of the Kleinman-Bylander projectors `Rₗₙ(|r|)Yₗₘ(r̂)` at a
given angular momentum `l`.
"""
function n_projector_radials(psp::AbstractPsP, l) end

"""
Number of radial parts of the pseudo-atomic wavefunctions with angular momentum `l`.
"""
function n_pseudo_orbital_radials(psp::AbstractPsP, l) end

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
function has_pseudo_orbitals(psp::AbstractPsP) end

"""
Projector coupling matrix for angular momentum `l`.
"""
function projector_coupling(psp::AbstractPsP, l) end

"""
Local part of the pseudopotential evaluated at real-space point `r`.
"""
function local_potential_real(psp::AbstractPsP)
     _ -> nothing
end

"""
The `n`th nonlocal Kleinman-Bylander projector at angular momentum `l` evaluated at
real-space point `r`.
"""
function projector_real(psp::AbstractPsP, l::Integer, n::Integer)
    _ -> nothing
end

"""
The `n`th pseudo-atomic orbital at angular momentum `l` evaulated at real-space point `r`.
"""
function pseudo_orbital_real(psp::AbstractPsP, l::Integer, n::Integer)
    _ -> nothing
end

"""
Pseudo-atomic valence charge density evaluated at real-space point `r`.
"""
function valence_charge_density_real(psp::AbstractPsP)
    _ -> nothing
end

"""
Model core charge density evaluated at real-space point `r`.
"""
function core_charge_density_real(psp::AbstractPsP)
    _ -> nothing
end

"""
Local part of the pseudopotential evaluated at reciprocal-space point `q`.
"""
function local_potential_fourier(psp::AbstractPsP)
    _ -> nothing
end

"""
The `n`th nonlocal Kleinman-Bylander projector at angular momentum `l` evaluated at
reciprocal-space point `q`.
"""
function projector_fourier(psp::AbstractPsP, l::Integer, n::Integer)
    _ -> nothing
end

"""
The `n`th pseudo-atomic orbital at angular momentum `l` evaulated at reciprocal-space
point `q`.
"""
function pseudo_orbital_fourier(psp::AbstractPsP, l::Integer, n::Integer)
    _ -> nothing
end

"""
Pseudo-atomic valence charge density evaluated at reciprocal-space point `q`.
"""
function valence_charge_density_fourier(psp::AbstractPsP)
    _ -> nothing 
end

"""
Model core charge density evaluated at reciprocal-space point `q`.
"""
function core_charge_density_fourier(psp::AbstractPsP)
    _ -> nothing
end

"""
Pseudopotential energy correction (the DC component of the Fourier transform of the
local part of the pseudopotential).
"""
function pseudo_energy_correction(T::Type, psp::AbstractPsP) end

#!!! Convenience functions !!!#
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
Number of radial parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_projector_radials(psp::AbstractPsP)
    return sum(l -> n_projector_radials(psp, l), 0:max_angular_momentum(psp); init=0)
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
    return sum(l -> n_projector_angulars(psp, l), 0:max_angular_momentum(psp); init=0)
end

"""
Number pseudo-atomic wavefunctions `Rₗₙ(|r|) * Yₗₘ(r̂)` at angular momenta `l` up to the
maximum angular momentum channel.
"""
function n_pseudo_orbital_radials(psp::AbstractPsP)
    return sum(l -> n_pseudo_orbital_radials(psp, l), 0:max_angular_momentum(psp); init=0)
end

"""
Number of angular parts of the pseudo-atomic wavefunctions with angular momentum `l`.
"""
function n_pseudo_orbital_angulars(psp::AbstractPsP, l)
    # for angular momentum l, magnetic q.n. m ∈ -l:l
    return n_pseudo_orbital_radials(psp::AbstractPsP, l) * (2l + 1)
end

"""
Number of angular parts of the pseudo-atomic wavefunctions at all angular momenta up
to the maximum angular momentum channel.
"""
function n_pseudo_orbital_angulars(psp::AbstractPsP)
    return sum(l -> n_pseudo_orbtial_angulars(psp, l), 0:max_angular_momentum(psp); init=0)
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

function Base.show(io::IO, psp::AbstractPsP)
    typename = string(typeof(psp))
    el = element(psp)
    z = valence_charge(psp)
    return print(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", psp::AbstractPsP)
    println(io, typeof(psp))
    @printf "%032s: %s\n" "formalism" formalism(psp)
    @printf "%032s: %s\n" "element" element(psp)
    @printf "%032s: %f\n" "valence charge" valence_charge(psp)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(psp)
    @printf "%032s: %s\n" "non-linear core correction" has_core_density(psp)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(psp)
    @printf "%032s: %s\n" "number of projectors" n_projector_radials(psp)
    @printf "%032s: %s" "number of pseudo-atomic orbitals" n_pseudo_orbital_radials(psp)
end
