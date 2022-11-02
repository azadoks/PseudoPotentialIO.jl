using DocStringExtensions

abstract type AbstractPsP end

"""
$(SIGNATURES)

Maximum angular momentum channel in the local part of the pseudopotential.
"""
function l_max(psp::AbstractPsP) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors `Rl(r)` at a given angular
momentum.
"""
function n_proj_radial(psp::AbstractPsP, l::Integer) end

"""
$(SIGNATURES)

Number of radial parts of the Kleinman-Bylander projectors `Rl(r)` at all angular momenta up
to the maximum angular momentum channel.
"""
n_proj_radial(psp::AbstractPsP) = sum(l -> n_proj_radial(psp, l), 0:l_max(psp); init=0)

"""
$(SIGNATURES)

Number of Kleinman-Bylander projectors `R(r) * Ylm(R)` at angular momentum `l`.
"""
n_proj(psp::AbstractPsP, l::Integer) = n_proj_radial(psp, l) * (2l + 1)

"""
$(SIGNATURES)

Number of Kleinman-Bylander projectors `R(r) * Ylm(R)` at angular momenta `l` up to the
maximum angular momentum channel.
"""
n_proj(psp::AbstractPsP) = n_proj_radial(psp) * (2l + 1)

"""
$(SIGNATURES)

Number of radial parts of the pseudo-atomic wavefunctions `χ(r)`.
"""
function n_pseudo_wfc(psp::AbstractPsP) end

"""
$(SIGNATURES)

Pseudo-atomic valence charge.
"""
function z_valence(psp::AbstractPsP) end

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

Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(psp::AbstractPsP) end

"""
$(SIGNATURES)

Whether the pseudopotential contains non-linear core correction data.
"""
function has_nlcc(psp::AbstractPsP) end

"""
$(SIGNATURES)

Kleinman-Bylander energy E_{ln,l'n'}.
"""
function e_kb(psp::AbstractPsP, li::Integer, ni::Integer, lj::Integer, nj::Integer) end

function e_kb(psp::AbstractPsP, l::Integer, n::Integer) end

"""
$(SIGNATURES)

Evaluate the local part of the pseudopotential at real-space radial distance `r` from the
nucleus.
"""
function v_local_real(psp::AbstractPsP, r) end

"""
$(SIGNATURES)

Evaluate the local part of the pseudopotential at real-space coordinate `R` in a coordinate
frame centered on the nucleus.
"""
v_local_real(psp::AbstractPsP, R::AbstractVector) = v_local_real(psp, norm(R))

"""
$(SIGNATURES)

Correction to the local part of the pseudopotential at real-space radial distance `r` from
the nucleus.
"""
v_local_correction_real(psp::AbstractPsP, r) = -z_valence(psp) / r

"""
$(SIGNATURES)

Fourier transform of the correction to the local part of the pseudopotential at
reciprocal-space distance `q` from the Γ-point.
"""
v_local_correction_fourier(psp::AbstractPsP, q) = 4π * -z_valence(psp) / q^2

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the local part of the pseudopotential at
the real-space radial distance `r` from the nucleus to the reciprocal space radial distance
`q` from the Γ-point.
"""
function v_local_ft_integrand(psp::AbstractPsP, q, r)
    return 4π * r^2 * sin(q * r) / (q * r) *
           (v_local_real(psp, r) - v_local_correction_real(psp, r))
end

"""
$(SIGNATURES)
"""
function v_local_fourier(psp::AbstractPsP, q)
    return quadgk(r -> v_local_ft_integrand(psp, q, r), 0, Inf)[1] +
           v_local_correction_fourier(psp, q)
end

function v_local_fourier(psp::AbstractPsP, G::AbstractVector)
    return v_local_fourier(psp, norm(G))
end

function projector_real(psp::AbstractPsP, l::Integer, n::Integer, r) end

function projector_real(psp::AbstractPsP, l::Integer, n::Integer,
                        R::AbstractVector)
    return projector_real(psp, l, n, norm(R))
end

function projector_ft_integrand(psp::AbstractPsP, l::Integer, n::Integer, q, r)
    return 4π * r^2 * sphericalbesselj(l, q * r) * projector_real(psp, l, n, r)
end

function projector_fourier(psp::AbstractPsP, l::Integer, n::Integer, q)
    return quadgk(r -> projector_ft_integrand(psp, l, n, q, r), 0, Inf)[1]
end

function projector_fourier(psp::AbstractPsP, l::Integer, n::Integer,
                           G::AbstractVector)
    return projector_fourier(psp, l, n, norm(G))
end

function _pseudo_energy_correction_integrand(psp::AbstractPsP, n_electrons::Integer, r)
    return 4π * n_electrons * r^2 * (v_local_real(psp, r) - v_local_correction_real(psp, r))
end

function pseudo_energy_correction(psp::AbstractPsP, n_electrons::Integer)
    return quadgk(r -> _pseudo_energy_correction_integrand(psp, n_electrons, r), 0, Inf)[1]
end
