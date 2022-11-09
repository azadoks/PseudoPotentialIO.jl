abstract type AbstractPsP end

Base.Broadcast.broadcastable(psp::AbstractPsP) = Ref(psp)

"""
$(SIGNATURES)

Element which the pseudopotential was constructed to reproduce.
"""
function element(psp::AbstractPsP)::PeriodicTable.Element end

"""
$(SIGNATURES)

Maximum angular momentum channel in the local part of the pseudopotential.
"""
function max_angular_momentum(psp::AbstractPsP) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors `Rl(r)` at a given angular
momentum.
"""
function n_projectors(psp::AbstractPsP, l) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors at all angular momenta up
to the maximum angular momentum channel.
"""
function n_projectors(psp::AbstractPsP)
    return sum(l -> n_projector_radials(psp, l), 0:max_angular_momentum(psp); init=0)
end

"""
$(SIGNATURES)

Number of radial parts of the pseudo-atomic wavefunctions with angular momentum `l`.
"""
function n_pseudo_orbitals(psp::AbstractPsP, l) end

"""
$(SIGNATURES)

Number pseudo-atomic wavefunctions `R(r) * Ylm(R)` at angular momenta `l` up to the maximum
angular momentum channel.
"""
function n_pseudo_orbitals(psp::AbstractPsP)
    return sum(n_pseudo_orbitals(psp, l), 0:max_angular_momentum(psp); init=0)
end

"""
$(SIGNATURES)

Pseudo-atomic valence charge.
"""
function valence_charge(psp::AbstractPsP) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(psp::AbstractPsP) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(psp::AbstractPsP) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(psp::AbstractPsP) end

"""
$(SIGNATURES)

Whether the pseudopotential is of the Coulombic kind.
"""
function is_coulomb(psp::AbstractPsP) end

"""
$(SIGNATURES)

Formalism of the pseudopotential (norm-conserving, ultrasoft, projector-augmented wave,
or Coulomb).
"""
function formalism(psp::AbstractPsP)::Symbol
    # The order here matters because some v2.0.1 PAW pseudos have both
    # is_paw and is_ultrasoft
    is_paw(psp) && return :paw
    is_ultrasoft(psp) && return :ultrasoft
    is_norm_conserving(psp) && return :norm_conserving
    is_coulomb(psp) && return :coulomb
end

"""
$(SIGNATURES)

Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(psp::AbstractPsP) end

"""
$(SIGNATURES)

Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(psp::AbstractPsP)::Symbol = has_spin_orbit(psp) ? :full : :scalar

"""
$(SIGNATURES)

Whether the pseudopotential contains non-linear core correction data (model core charge
density).
"""
function has_nlcc(psp::AbstractPsP) end

"""
$(SIGNATURES)

Projector coupling constant between the `ni`th projector with angular momentum `li` and the
`nj`th projector with angular momentum `lj`. Sometimes called the Kleinman-Bylander energy,
depending on the choice of units.
"""
function projector_coupling(psp::AbstractPsP, li, ni, lj, nj) end

"""
$(SIGNATURES)

Projector coupling constant between the `n`th projector with angular momentum `l` itself.
Sometimes called the Kleinman-Bylander energy, depending on the choice of units.
"""
function projector_coupling(psp::AbstractPsP, l, n) end

"""
$(SIGNATURES)

Local part of the pseudopotential at real-space point `r`.
"""
function local_potential_real(psp::AbstractPsP, r::Real) end

"""
$(SIGNATURES)

Correction to the local part of the pseudopotential at real-space point `r`.
"""
function local_potential_correction_real(psp::AbstractPsP, r::T)::T where {T<:Real}
    return -valence_charge(psp) / r
end

"""
$(SIGNATURES)

Fourier transform of the correction to the local part of the pseudopotential at
reciprocal-space point `q`.
"""
function local_potential_correction_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return 4T(π) * -valence_charge(psp) / q^2
end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the local part of the pseudopotential at
the real-space point `r` to the reciprocal space point `q`.
"""
function local_potential_fourier_transform_integrand(psp::AbstractPsP, q::T,
                                                     r::T)::T where {T<:Real}
    Vlocal_corrected = (local_potential_real(psp, r) -
                        local_potential_correction_real(psp, r))
    return 4T(π) * r^2 * sphericalbesselj(0, q * r) * Vlocal_corrected
end

"""
$(SIGNATURES)

Fourier transform of the local part of the pseudopotential at reciprocal-space point `q`.
"""
function local_potential_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return quadgk(r -> local_potential_fourier_transform_integrand(psp, q, r), 0, Inf)[1] +
           local_potential_correction_fourier(psp, q)
end

"""
$(SIGNATURES)

Value of the `n`th Kleinman-Bylander non-local projector with angular momentum `l` at
real-space point `r`.
"""
function projector_radial_real(psp::AbstractPsP, l, n, r) end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the `n`th Kleinman-Bylander non-local
projector with angular momentum `l` at real-space point `r` and reciprocal-space point `q`.
"""
function projector_radial_fourier_transform_integrand(psp::AbstractPsP, l::Int, n::Int,
                                                      q::T,
                                                      r::T)::T where {T<:Real}
    return 4π * r^2 * sphericalbesselj(l, q * r) * projector_radial_real(psp, l, n, r)
end

"""
$(SIGNATURES)

Fourier transform of the `n`th Kleinman-Bylander non-local projector with angular momentum
`l` at reciprocal-space point `q`.
"""
function projector_radial_fourier(psp::AbstractPsP, l, n, q::T)::T where {T<:Real}
    return quadgk(r -> projector_radial_fourier_transform_integrand(psp, l, n, q, r), 0,
                  Inf)[1]
end

"""
$(SIGNATURES)

Integrand for the "DC" component of the local part of the pseudopotential.
"""
function pseudo_energy_correction_fourier_transform_integrand(psp::AbstractPsP,
                                                              r::T)::T where {T<:Real}
    return 4T(π) * r^2 *
           (local_potential_real(psp, r) - local_potential_correction_real(psp, r))
end

"""
$(SIGNATURES)

"DC" component of the local part of the pseudopotential.
"""
function pseudo_energy_correction(psp::AbstractPsP)
    return quadgk(r -> pseudo_energy_correction_fourier_transform_integrand(psp, r), 0,
                  Inf)[1]
end

"""
$(SIGNATURES)

Model core charge density at real-space point `r`.
"""
function core_charge_density_real(psp::AbstractPsP, r::Real) end

"""
$(SIGNATURES)

Model core charge density at real-space vector `r`.
"""
function core_charge_density_real(psp::AbstractPsP, R::AbstractVector{T})::T where {T<:Real}
    return core_charge_density_real(psp, norm(R))
end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the model core charge density at
the real-space point `r` to the reciprocal space point `q`.
"""
function core_charge_density_fourier_transform_integrand(psp::AbstractPsP, q::T,
                                                         r::T)::T where {T<:Real}
    return 4T(π) * r^2 * sphericalbesselj_fast(0, q * r) * core_charge_density_real(psp, r)
end

"""
$(SIGNATURES)

Fourier transform of the model core charge density at reciprocal-space point `q`.
"""
function core_charge_density_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return quadgk(r -> core_charge_density_fourier_transform_integrand(psp, q, r), 0, Inf)[1]
end

"""
$(SIGNATURES)

Fourier transform of the model core charge density at reciprocal-space vector `K`.
"""
function core_charge_density_fourier(psp::AbstractPsP,
                                     K::AbstractVector{T})::T where {T<:Real}
    return core_charge_density_fourier(psp, norm(K))
end

"""
$(SIGNATURES)

Valence charge density at real-space point `r`.
"""
function valence_charge_density_real(psp::AbstractPsP, r::Real) end

"""
$(SIGNATURES)

Valence charge density at real-space vector `r`.
"""
function valence_charge_density_real(psp::AbstractPsP,
                                     R::AbstractVector{T})::T where {T<:Real}
    return valence_charge_density_real(psp, norm(R))
end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the valence charge density at
the real-space point `r` to the reciprocal space point `q`.
"""
function valence_charge_density_fourier_transform_integrand(psp::AbstractPsP, q::T,
                                                            r::T)::T where {T<:Real}
    return 4T(π) * r^2 * sphericalbesselj_fast(0, q * r) *
           valence_charge_density_real(psp, r)
end

"""
$(SIGNATURES)

Fourier transform of the valence charge density at reciprocal-space point `q`.
"""
function valence_charge_density_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return quadgk(r -> valence_charge_density_fourier_transform_integrand(psp, q, r), 0,
                  Inf)[1]
end

"""
$(SIGNATURES)

Fourier transform of the valence charge density at reciprocal-space vector `K`.
"""
function valence_charge_density_fourier(psp::AbstractPsP,
                                        K::AbstractVector{T})::T where {T<:Real}
    return valence_charge_density_fourier(psp, norm(K))
end

function Base.show(io::IO, psp::AbstractPsP)
    typename = string(typeof(psp))
    el = element(psp).symbol
    z = valence_charge(psp)
    return println(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", psp::AbstractPsP)
    println(io, typeof(psp))
    @printf "%032s: %s\n" "formalism" formalism(psp)
    @printf "%032s: %s\n" "element" element(psp).symbol
    @printf "%032s: %f\n" "valence charge" valence_charge(psp)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(psp)
    @printf "%032s: %s\n" "non-linear core correction" has_nlcc(psp)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(psp)
    @printf "%032s: %s\n" "number of projectors" n_projectors(psp)
    @printf "%032s: %s"   "number of pseudo-atomic orbitals" n_pseudo_orbitals(psp)
end
