abstract type AbstractPsP end

"""
$(SIGNATURES)

Pseudopotential file format.
"""
function format(psp::AbstractPsP) end

"""
$(SIGNATURES)

Element which the pseudopotential was constructed to reproduce.
"""
function element(psp::AbstractPsP)::PeriodicTable.Element end

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
function n_proj_radial(psp::AbstractPsP, l) end

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
n_proj(psp::AbstractPsP, l) = n_proj_radial(psp, l) * (2l + 1)

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

Formalism of the pseudopotential (norm-conserving, ultrasoft, projector-augmented wave,
or Coulomb).
"""
function formalism(psp::AbstractPsP)::Symbol
    is_norm_conserving(psp) && return :norm_conserving
    is_ultrasoft(psp) && return :ultrasoft
    is_paw(psp) && return :projector_augmented_wave
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

Whether the pseudopotential contains non-linear core correction data.
"""
function has_nlcc(psp::AbstractPsP) end

"""
$(SIGNATURES)

Kleinman-Bylander energy E_{ln,l'n'}.
"""
function e_kb(psp::AbstractPsP, li, ni, lj, nj) end

function e_kb(psp::AbstractPsP, l, n) end

"""
$(SIGNATURES)

Evaluate the local part of the pseudopotential at real-space point `r`.
"""
function v_local_real(psp::AbstractPsP, r::Real) end

"""
$(SIGNATURES)

Evaluate the local part of the pseudopotential at real-space vector `R`.
"""
function v_local_real(psp::AbstractPsP, R::AbstractVector{T})::T where {T<:Real}
    return v_local_real(psp, norm(R))
end

"""
$(SIGNATURES)

Correction to the local part of the pseudopotential at real-space point `r`.
"""
function v_local_correction_real(psp::AbstractPsP, r::T)::T where {T<:Real}
    return -z_valence(psp) / r
end

"""
$(SIGNATURES)

Fourier transform of the correction to the local part of the pseudopotential at
reciprocal-space point `q`.
"""
function v_local_correction_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return 4T(π) * -z_valence(psp) / q^2
end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the local part of the pseudopotential at
the real-space point `r` to the reciprocal space point `q`.
"""
function v_local_ft_integrand(psp::AbstractPsP, q::T, r::T)::T where {T<:Real}
    v_local_corrected = (v_local_real(psp, r) - v_local_correction_real(psp, r))
    return 4T(π) * r^2 * sphericalbesselj(0, q * r) * v_local_corrected
end

"""
$(SIGNATURES)

Fourier transform of the local part of the pseudopotential at reciprocal-space point `q`.
"""
function v_local_fourier(psp::AbstractPsP, q::T)::T where {T<:Real}
    return quadgk(r -> v_local_ft_integrand(psp, q, r), 0, Inf)[1] +
           v_local_correction_fourier(psp, q)
end

"""
$(SIGNATURES)

Fourier transform of the local part of the pseudopotential at reciprocal-space vector `K`.
"""
function v_local_fourier(psp::AbstractPsP, K::AbstractVector{T})::T where {T<:Real}
    return v_local_fourier(psp, norm(K))
end

"""
$(SIGNATURES)

Value of the `n`th Kleinman-Bylander non-local projector with angular momentum `l` at
real-space point `r`.
"""
function projector_real(psp::AbstractPsP, l, n, r) end

"""
$(SIGNATURES)

Value of the `n`th Kleinman-Bylander non-local projector with angular momentum `l` at
real-space vector `R`.
"""
function projector_real(psp::AbstractPsP, l::Int, n::Int,
                        R::AbstractVector{T})::T where {T<:Real}
    return projector_real(psp, l, n, norm(R))
end

"""
$(SIGNATURES)

Integrand for the centered Fourier transform of the `n`th Kleinman-Bylander non-local
projector with angular momentum `l` at real-space point `r` and reciprocal-space point `q`.
"""
function projector_ft_integrand(psp::AbstractPsP, l::Int, n::Int, q::T,
                                r::T)::T where {T<:Real}
    return 4π * r^2 * sphericalbesselj(l, q * r) * projector_real(psp, l, n, r)
end

"""
$(SIGNATURES)

Fourier transform of the `n`th Kleinman-Bylander non-local projector with angular momentum
`l` at reciprocal-space point `q`.
"""
function projector_fourier(psp::AbstractPsP, l, n, q)
    return quadgk(r -> projector_ft_integrand(psp, l, n, q, r), 0, Inf)[1]
end

"""
$(SIGNATURES)
"""
function projector_fourier(psp::AbstractPsP, l, n,
                           G::AbstractVector{T}) where {T<:Real}
    return projector_fourier(psp, l, n, norm(G))
end

"""
$(SIGNATURES)
"""
function pseudo_energy_correction_ft_integrand(psp::AbstractPsP, r)
    return 4π * r^2 * (v_local_real(psp, r) - v_local_correction_real(psp, r))
end

"""
$(SIGNATURES)
"""
function pseudo_energy_correction(psp::AbstractPsP)
    return quadgk(r -> pseudo_energy_correction_ft_integrand(psp, r), 0, Inf)[1]
end

function Base.show(io::IO, psp::AbstractPsP)
    typename = string(typeof(psp))
    el = element(psp).symbol
    z = z_valence(psp)
    println(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", psp::AbstractPsP)
    println(io, format(psp))
    @printf "%032s: %s\n" "formalism" formalism(psp)
    @printf "%032s: %s\n" "element" element(psp).symbol
    @printf "%032s: %f\n" "valence charge" z_valence(psp)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(psp)
    @printf "%032s: %s\n" "non-linear core correction" has_nlcc(psp)
    @printf "%032s: %d\n" "maximum angular momentum" l_max(psp)
    @printf "%032s: %s\n" "number of projectors" n_proj_radial(psp)
    @printf "%032s: %s" "number of pseudo-atomic orbitals" n_pseudo_wfc(psp)
end
