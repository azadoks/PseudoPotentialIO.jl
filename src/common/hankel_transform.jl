@doc raw"""
Hankel / Bessel-Fourier transform.
The function to transform should be rapidly decaying to zero within the bounds of the mesh.
Note that the `hankel_transform` requires `r²f` as input, i.e. the function must be
pre-multiplied by the square of the radial mesh.

The transform is defined as:

```math
4\pi \int_0^{\infty} r^2 f(r) j_l(q r) dr \approx
4\pi \int_{r_1}^{r_N} r^2 f(r) j_l(q r) dr
```

where ``j_l(x)`` is the spherical Bessel function of the first kind at order ``l``.
"""
function hankel_transform end

@inbounds function hankel_transform(mesh::RadialMesh{T}, r²f::AbstractVector{T},
                                    qs::AbstractVector{S}, l::Integer,
                                    weights::AbstractVector{S},
                                    integrand::AbstractVector{S};
                                    quadrature_method=Simpson()) where {T<:Real,S<:Real}
    r = @view mesh[1:length(r²f)]
    weights = @view weights[1:length(r²f)]
    integrand = @view integrand[1:length(r²f)]

    integration_weights!(weights, mesh, quadrature_method)
    jₗ = fast_sphericalbesselj(l)

    f = map(qs) do q
        integrand .= r²f .* jₗ.(q .* r)
        integral = dot(weights, integrand)
        return 4π * integral
    end

    return f
end

function hankel_transform(psp::NormConservingPsP{T}, qs::AbstractVector{S};
                          quadrature_method=Simpson(),
                          correction_method=CoulombCorrection())::NormConservingPsP{S} where {T<:Real,
                                                                                              S<:Real}
    weights = Vector{S}(undef, length(psp.r))
    integrand = Vector{S}(undef, length(psp.r))

    Vloc = hankel_transform(psp, qs, weights, integrand, NumericLocalPotential();
                            quadrature_method, correction_method)
    β = hankel_transform(psp, qs, weights, integrand, NumericProjector(); quadrature_method)
    χ = hankel_transform(psp, qs, weights, integrand, NumericState(); quadrature_method)
    ρval = hankel_transform(psp, qs, weights, integrand, ValenceChargeDensity();
                            quadrature_method)
    ρcore = hankel_transform(psp, qs, weights, integrand, CoreChargeDensity();
                             quadrature_method)

    return NormConservingPsP{S}(psp.identifier, nothing, psp.Zatom, psp.Zval, psp.lmax, qs,
                             Vloc, β, psp.D, χ, ρcore, ρval)
end

function hankel_transform(psp::UltrasoftPsP{T}, qs::AbstractVector{S};
                          quadrature_method=Simpson(),
                          correction_method=CoulombCorrection())::UltrasoftPsP{S} where {T<:Real,
                                                                                         S<:Real}
    weights = Vector{S}(undef, length(psp.r))
    integrand = Vector{S}(undef, length(psp.r))

    Vloc = hankel_transform(psp, qs, weights, integrand, NumericLocalPotential();
                            quadrature_method, correction_method)
    β = hankel_transform(psp, qs, weights, integrand, NumericProjector(); quadrature_method)
    χ = hankel_transform(psp, qs, weights, integrand, NumericState(); quadrature_method)
    Q = hankel_transform(psp, qs, weights, integrand, AugmentationFunction();
                         quadrature_method)
    ρval = hankel_transform(psp, qs, weights, integrand, ValenceChargeDensity();
                            quadrature_method)
    ρcore = hankel_transform(psp, qs, weights, integrand, CoreChargeDensity();
                             quadrature_method)

    # TODO: broken, type of Χ gives a MethodError
    return UltrasoftPsP{S}(psp.identifier, nothing, psp.Zatom, psp.Zval, psp.lmax, qs, Vloc, β,
                           psp.D, χ, Q, psp.q, ρcore, ρval)
end

#TODO: remove the duplication here by having get_quantity take args...
function hankel_transform(psp::NumericPsP{T}, qs::AbstractVector{S},
                          weights::AbstractVector{S},
                          integrand::AbstractVector{S}, quantity::AbstractPsPQuantity;
                          quadrature_method=Simpson())::Union{Nothing,
                                                              Vector{S}} where {T<:Real,
                                                                                S<:Real}
    !has_quantity(psp, quantity) && return nothing
    r²f = get_quantity(psp, quantity)
    return hankel_transform(psp.r, r²f, qs, 0, weights, integrand; quadrature_method)
end

function hankel_transform(psp::NumericPsP{T}, qs::AbstractVector{S},
                          weights::AbstractVector{S},
                          integrand::AbstractVector{S}, quantity::AbstractPsPQuantity, l::Integer,
                          args...;
                          quadrature_method=Simpson())::Union{Nothing,
                                                              Vector{S}} where {T<:Real,
                                                                                S<:Real}
    !has_quantity(psp, quantity) && return nothing
    r²f = get_quantity(psp, quantity, l, args...)
    return hankel_transform(psp.r, r²f, qs, l, weights, integrand; quadrature_method)
end

function hankel_transform(psp::NumericPsP{T}, qs::AbstractVector{S},
                          weights::AbstractVector{S},
                          integrand::AbstractVector{S}, quantity::NumericLocalPotential;
                          quadrature_method=Simpson(),
                          correction_method=CoulombCorrection())::Vector{S} where {T<:Real,
                                                                                   S<:Real}
    Vloc = get_quantity(psp, quantity)

    r = @view psp.r[1:length(Vloc)]
    weights = @view weights[1:length(Vloc)]
    integrand = @view integrand[1:length(Vloc)]

    rVloc_m_corr = r .* Vloc .-
                   local_potential_correction.(Ref(psp), Ref(correction_method),
                                               Ref(RealSpace()), r)

    integration_weights!(weights, psp.r, quadrature_method)

    f = map(qs) do q
        integrand .= rVloc_m_corr .* sin.(q .* r)
        integral = dot(weights, integrand)
        corr = local_potential_correction(psp, correction_method, FourierSpace(), q)
        return 4π * (integral / q + corr)
    end

    return f
end

function hankel_transform(psp::NumericPsP{T}, qs::AbstractVector{S},
                          weights::AbstractVector{S},
                          integrand::AbstractVector{S}, quantity::AbstractProjector;
                          quadrature_method=Simpson()) where {T<:Real,S<:Real}
    if has_quantity(psp, quantity)
        f = map(angular_momenta(psp)) do l
            map(1:n_radials(psp, quantity, l)) do n
                return hankel_transform(psp, qs, weights, integrand, quantity, l, n;
                                        quadrature_method)
            end
        end
        f = OffsetArray(f, angular_momenta(psp))
    else
        f = nothing
    end

    return f
end

# TODO: this is most likely wrong!
function hankel_transform(psp::UltrasoftPsP{T}, qs::AbstractVector{S},
                          weights::AbstractVector{S},
                          integrand::AbstractVector{S}, quantity::AugmentationFunction;
                          quadrature_method=Simpson()) where {T<:Real,S<:Real}
    Q = map(angular_momenta(psp)) do l
        Ql = similar(psp.Q[l])
        for n in axes(psp.Q[l], 1)
            for m in axes(psp.Q[l], 2)
                Ql[n, m] = hankel_transform(psp, qs, weights, integrand, quantity, l, n, m;
                                            quadrature_method)
            end
        end
        return Ql
    end
    Q = OffsetArray(Q, angular_momenta(psp))

    return Q
end
