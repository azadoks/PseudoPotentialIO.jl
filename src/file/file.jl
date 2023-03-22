@doc raw"""
Abstract type representing a pseudopotential file.

The structure of the data should closely mirror the format of the file, and the values
of quantities should be exactly those found in the file.

Required methods:

```julia
# A unique string, usually a hash or checksum.
function identifier(file::PsPFile)::AbstractString end
# A short string listing the file format (e.g. `"PSP8"`)
function format(file::PsPFile)::AbstractString end
# The symbol of the element for which the file contains a pseudopotential (e.g. `"Ag"`)
function elemental_symbol(file::PsPFile)::AbstractString end
# The maximum angular momentum channel contained in the file
function max_angular_momentum(file::PsPFile)::Integer end
# The number of non-local projectors for angular momentum `l` contained in the file
function n_beta_projector_radials(file::PsPFile, l::Integer)::Integer end
# The number of chi functions for angular momentum `l` contained in the file
function n_chi_projector_radials(file::PsPFile, l::Integer)::Integer end
# The pseudo-atomic valence charge
function valence_charge(file::PsPFile)::Real end
# Whether the file contains a norm-conserving pseudopotential
function is_norm_conserving(file::PsPFile)::Bool end
# Whether the file contains an ultrasoft pseudopotential
function is_ultrasoft(file::PsPFile)::Bool end
# Whether the file contains a projector-augmented wave pseudopotential
function is_paw(file::PsPFile)::Bool end
# Whether the file contains a pseudopotential supporting spin-orbit coupled calculations
function has_spin_orbit(file::PsPFile)::Bool end
# Whether the file contains a pseudopotential supporting non-linear core corrections
function has_core_density(file::PsPFile)::Bool end
```
"""
abstract type PsPFile end

#!!! Required functions !!!#
"""
Identifying data (preferably unique).
"""
function identifier(file::PsPFile) end

"""
Pseudopotential file format.
"""
function format(file::PsPFile) end

"""
The element which the pseudopotential was constructed to reproduce.
"""
function element(file::PsPFile) end

"""
Maximum angular momentum channel in the local part of the pseudopotential.
"""
function max_angular_momentum(file::PsPFile) end

"""
Number of radial functions of the given quantity contained in the pseudopotential file.
"""
function n_radials(::PsPProjector, file::PsPFile, l) end

"""
Pseudo-atomic valence charge.
"""
function valence_charge(file::PsPFile) end

"""
Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(file::PsPFile) end

"""
Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(file::PsPFile) end

"""
Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(file::PsPFile) end

"""
Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(file::PsPFile) end

#!!! Convenience functions !!!#
"""
Formalism of the pseudopotential.
"""
function formalism(file::PsPFile)::Type
    is_paw(file) && return ProjectorAugmentedWavePsP
    is_ultrasoft(file) && return UltrasoftPsP
    is_norm_conserving(file) && return NormConservingPsP
end

"""
Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(file::PsPFile)::Symbol = has_spin_orbit(file) ? :full : :scalar

function n_radials(q::PsPProjector, file::PsPFile)
    return sum(l -> n_radials(q, file, l), 0:max_angular_momentum(file); init=0)
end

function n_angulars(q::PsPProjector, file::PsPFile, l)
    return n_radials(q, file, l) * (2l + 1)
end

"""
Number of angular parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_angulars(q::PsPProjector, file::PsPFile)
    return sum(l -> n_angulars(q, file, l), 0:max_angular_momentum(file); init=0)
end

Base.Broadcast.broadcastable(file::PsPFile) = Ref(file)

Base.:(==)(file1::PsPFile, file2::PsPFile) = file1.checksum == file2.checksum
Base.hash(file::PsPFile) = hash(file.checksum)

function Base.show(io::IO, file::PsPFile)
    typename = string(typeof(file))
    el = element(file)
    z = valence_charge(file)
    return print(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", file::PsPFile)
    println(io, typeof(file))
    @printf "%032s: %s\n" "identifier" identifier(file)
    @printf "%032s: %s\n" "formalism" formalism(file)
    @printf "%032s: %s\n" "element" element(file)
    @printf "%032s: %f\n" "valence charge" valence_charge(file)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(file)
    @printf "%032s: %s\n" "non-linear core correction" has_quantity(CoreChargeDensity(), file)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(file)
    @printf "%032s: %s\n" "number of beta projectors" n_radials(BetaProjector(), file)
    @printf "%032s: %s" "number of chi projectors" n_radials(ChiProjector(), file)
end
