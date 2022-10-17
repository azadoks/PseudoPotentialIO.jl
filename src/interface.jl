abstract type AbstractPsP end

"""
Maximum angular momentum channel in the local part of the pseudopotential.
"""
function ang_mom_max(psp::AbstractPsP) end

"""
Angular momentum channel of the local part of the pseudopotential.
"""
function l_local(psp::AbstractPsP) end

"""
Number of radial parts of the Kleinman-Bylander projectors `R(r)` at a given angular
momentum.
"""
function n_proj_radial(psp::AbstractPsP, ang_mom::Integer) end

"""
Number of radial parts of the Kleinman-Bylander projectors `R(r)` at all angular momenta up
to the maximum angular momentum channel.
"""
function n_proj_radial(psp::AbstractPsP)
	return sum(l -> n_proj_radial(psp, l), 0:ang_mom_max(psp); init=0)
end

"""
Number of Kleinman-Bylander projectors `R(r) * Ylm(R)` at angular momentum `l`.
"""
function n_proj(psp::AbstractPsp, ang_mom::Integer)
	return n_proj_radial(psp, ang_mom) * (2 * ang_mom + 1)
end

"""
Number of Kleinman-Bylander projectors `R(r) * Ylm(R)` at angular momenta `l` up to the
maximum angular momentum channel.
"""
function n_proj(psp::AbstractPsp)
	return n_proj_radial(psp) * (2 * ang_mom + 1)
end

"""
Number of radial parts of the pseudo-atomic wavefunctions `χ(r)`.
"""
function n_wfc(psp::AbstractPsP) end

"""
Pseudo-atomic valence charge.
"""
function z_valence(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(psp::AbstractPsP) end

"""
Whether the pseudopotential is of the Coulombic kind.
"""
function is_coulomb(psp::AbstractPsP) end

"""
Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(psp::AbstractPsP) end

"""
Whether the pseudopotential contains non-linear core correction data.
"""
function has_nlcc(psp::AbstractPsP) end

"""
Evaluate the local part of the pseudopotential at real-space radial distance `r` from the
nucleus.
"""
function v_local_real(psp::AbstractPsP, r) end

"""
Evaluate the local part of the pseudopotential at real-space coordinate `R` in a coordinate
frame centered on the nucleus.
"""
function v_local_real(psp::AbstractPsP, R::AbstractVector)
	return v_local_real(psp, norm(R))
end

function v_local_correction_real(psp::AbstractPsP, r)
	return -z_valence(psp) / r
end

function v_local_correction_fourier(psp::AbstractPsP, q)
	return 4π * -z_valence(psp) / q^2
end

"""
Integrand for the centered Fourier transform of the local part of the pseudopotential at
the real-space radial distance `r` from the nucleus to the reciprocal space radial distance
`q` from the Γ-point.

4π r^2 sin(qr) / (qr) (V_{local}(r))

This expression is not well-localized because the local potential contains a Coulomb-like
term. To localize the integrand, a Coulomb-like correction is subtracted from the local part
of the pseudopotential:

4π r^2 sin(qr) / (qr) (V_{local}(r) - C_{local}(r)).
"""
function _v_local_ft_integrand(psp::AbstractPsp, q, r)
	return 4π * r^2 * sin(q * r) / (q * r) *
		   (v_local_real(psp, r) - v_local_correction_real(psp, r))
end

"""
Evaluate the local part of the pseudopotential at reciprocal-space distance `q` from the
Γ-point.

V_{local}(q) = ∫_{R+} 4π r^2 sin(qr) / (qr) (V_{local}(r) - C_{local}(r)) + F[C_{local}]
"""
function v_local_fourier(psp::AbstractPsP, q)
	return quadgk(r -> _v_local_ft_integrand(psp, q, r), 0, Inf)[1] +
		   v_local_correction_fourier(psp, q)
end

function v_local_fourier(psp::AbstractPsP, G::AbstractVector)
	return v_local_fourier(psp, norm(G))
end

function projector_real(psp::AbstractPsP, ang_mom::Integer, n::Integer, r) end

function projector_real(psp::AbstractPsP, ang_mom::Integer, n::Integer,
						R::AbstractVector)
	return projector_real(psp, ang_mom, n, norm(R))
end

function _projector_ft_integrand(psp::AbstractPsP, ang_mom::Integer, n::Integer, q, r)
	return 4π * r^2 * sphericalbesselj(ang_mom, q * r) * projector_real(psp, ang_mom, n, r)
end

function projector_fourier(psp::AbstractPsP, ang_mom::Integer, n::Integer, q)
	return quadgk(r -> _projector_ft_integrand(psp, ang_mom, n, q, r), 0, Inf)[1]
end

function projector_fourier(psp::AbstractPsP, ang_mom::Integer, n::Integer,
						   G::AbstractVector)
	return projector_fourier(psp, ang_mom, n, norm(G))
end

function _pseudo_energy_correction_integrand(psp::AbstractPsP, n_electrons::Integer, r)
	return 4π * n_electrons * r^2 * (v_local_real(psp, r) + z_valence(psp) / r)
end

function pseudo_energy_correction(psp::AbstractPsP, n_electrons::Integer)
	return quadgk(r -> _pseudo_energy_correction_integrand(psp, n_electrons, r), 0, Inf)[1]
end
