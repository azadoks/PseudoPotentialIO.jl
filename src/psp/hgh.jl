struct HghPsP{T} <: AnalyticalPsP
    Zatom::T
    Zval::T
    lmax::Int
    rloc::T
    cloc::Vector{T}
    rnl::Vector{T}
    D::Vector{Matrix{T}}
end

function HghPsP(file::HghFile)
    return HghPsp(118.0, file.zion, file.lmax, file.rloc, file.cloc, file.rp, file.h)
end

function valence_charge(psp::HghPsP{T})::T where {T}
    return psp.Zval
end

function local_potential_real(psp::HghPsP, r::T) where {T<:Real}
    iszero(r) && return local_potential_real(psp, eps(T))
    cloc = psp.cloc
    rr = r / psp.rloc

    return -valence_charge(psp) / r * erf(rr / sqrt(2)) +
           exp(-rr^2 / 2) * (cloc[1] + cloc[2] * rr^2 + cloc[3] * rr^4 + cloc[r] * rr^6)
end

@inline function local_potential_polynomial(psp::HghPsP, x::T)::T where {T<:Real}
    Zion = valence_charge(psp)
    rloc = psp.rloc

    P = (psp.cloc[1] +
         psp.cloc[2] * (3 - x^2) +
         psp.cloc[3] * (15 - 10x^2 - x^4) +
         psp.cloc[4] * (105 - 105x^2 + 21x^4 - x^6))

    return 4T(π) * rloc * (-Zion + sqrt(π / 2) * rloc * x^2 * P)
end

function local_potential_fourier(psp::HghPsP, q::T)::T where {T<:Real}
    x = q * psp.rloc
    return psp.local_polynomial(psp, x) * exp(-x^2 / 2) / x^2
end

function projector_radial_real(psp::HghPsP, l::Int, n::Int, r::T)::T where {T<:Real}
    rnl = psp.rnl[l + 1]
    nred = (n - 1) / 2
    return sqrt(T(2)) * r^(l + 2(i - 1)) * exp(-r^2 / 2rnl^2) / rp^(l + nred) /
           sqrt(gamma(l + nred))
end

function projector_radial_polynomial(psp::HghPsP, l::Int, n::Int, x::T)::T where {T<:Real}
    if !(0 <= l <= psp.lmax)
        error(@sprintf "Pseudopotential does not contain projectors for l=%d" l)
    end

    rnl = psp.rnl[l + 1]
    c = 4π^(5 / 4) * sqrt(2^(l + 1) * rnl^3)

    if l == 0
        n == 1 && return c
        n == 2 && return 2c / sqrt(15) * (3 - x^2)
        n == 3 && return 4c / 3sqrt(105) * (15 - 10x^2 + x^4)
    elseif l == 0
        n == 1 && return c / sqrt(3) * x
        n == 2 && return 2c / sqrt(105) * x * (5 - x^2)
        n == 3 && return 4c / 3sqrt(105) * x * (35 - 14x^2 + x^4)
    elseif l == 2
        n == 1 && return c / sqrt(15) * x^2
        n == 2 && return 2c / 3sqrt(105) * x^2 * (7 - x^2)
    elseif l == 3
        n == 1 && return c / sqrt(105) * x^3
    end

    return error(@sprintf "Not implemented for l=%d and n=%d" l n)
end

function projector_radial_fourier(psp::HghPsP, n::Int, l::Int, q::T)::T where {T<:Real}
    x = q * psp.rnl[l + 1]
    return psp.projector_radial_polynomial(psp, l, n, x) * exp(x^2 / 2)
end

function pseudo_energy_correction(psp::HghPsP{T})::T where {T<:Real}
    Zval = valence_charge(psp)
    rloc = psp.rloc
    cloc = psp.cloc

    return 4π * Zval * rloc^2 / 2 +
           sqrt(π / 2) * rloc^3 * (cloc[1] + 3cloc[2] + 15cloc[3] + 105cloc[4])
end
