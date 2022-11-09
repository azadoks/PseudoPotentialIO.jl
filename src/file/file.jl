abstract type PsPFile end

Base.Broadcast.broadcastable(file::PsPFile) = Ref(file)

"""
$(SIGNATURES)

Pseudopotential file format.
"""
function format(file::PsPFile) end

"""
$(SIGNATURES)

Element which the pseudopotential was constructed to reproduce.
"""
function element(file::PsPFile)::PeriodicTable.Element end

"""
$(SIGNATURES)

Maximum angular momentum channel in the local part of the pseudopotential.
"""
function max_angular_momentum(file::PsPFile) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors `Rl(r)` at a given angular
momentum.
"""
function n_projectors(file::PsPFile, l) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_projectors(file::PsPFile)
    return sum(l -> n_projector_radials(file, l), 0:max_angular_momentum(file); init=0)
end

"""
$(SIGNATURES)

Number of radial parts of the pseudo-atomic wavefunctions with angular momentum `l`.
"""
function n_pseudo_orbitals(file::PsPFile, l) end

"""
$(SIGNATURES)

Number pseudo-atomic wavefunctions `R(r) * Ylm(R)` at angular momenta `l` up to the maximum
angular momentum channel.
"""
function n_pseudo_orbitals(file::PsPFile)
    return sum(n_pseudo_orbitals(file, l), 0:max_angular_momentum(file); init=0)
end

"""
$(SIGNATURES)

Pseudo-atomic valence charge.
"""
function valence_charge(file::PsPFile) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(file::PsPFile) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(file::PsPFile) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(file::PsPFile) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the Coulombic kind.
"""
function is_coulomb(file::PsPFile) end

"""
$(SIGNATURES)

Formalism of the pseudopotential (norm-conserving, ultrasoft, projector-augmented wave,
or Coulomb).
"""
function formalism(file::PsPFile)::Symbol
    # The order here matters because some v2.0.1 PAW pseudos have both
    # is_paw and is_ultrasoft
    is_paw(file) && return :paw
    is_ultrasoft(file) && return :ultrasoft
    is_norm_conserving(file) && return :norm_conserving
    is_coulomb(file) && return :coulomb
end

"""
$(SIGNATURES)

Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(file::PsPFile) end

"""
$(SIGNATURES)

Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(file::PsPFile)::Symbol = has_spin_orbit(file) ? :full : :scalar

"""
$(SIGNATURES)

Whether the pseudopotential contains non-linear core correction data (model core charge
density).
"""
function has_nlcc(file::PsPFile) end

function Base.show(io::IO, file::PsPFile)
    typename = string(typeof(file))
    el = element(file).symbol
    z = valence_charge(file)
    return println(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", file::PsPFile)
    println(io, typeof(file))
    @printf "%032s: %s\n" "formalism" formalism(file)
    @printf "%032s: %s\n" "element" element(file).symbol
    @printf "%032s: %f\n" "valence charge" valence_charge(file)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(file)
    @printf "%032s: %s\n" "non-linear core correction" has_nlcc(file)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(file)
    @printf "%032s: %d\n" "number of projectors" n_projectors(file)
    @printf "%032s: %d"   "number of pseudo-atomic orbitals" n_pseudo_orbitals(file)
end