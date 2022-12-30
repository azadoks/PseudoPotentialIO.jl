#TODO Implement in-place versions which take a working array to use for storing
#TODO the integrand
@doc raw"""
Bessel-Fourier (also Hankel) transform of order `l` of a function on a radial mesh.
The function should be rapidly decaying to zero within the bounds of the mesh.

```math
4\pi \int_0^\inf f(r) j_l(q r) r^2 dr \approx
4\pi \int_{r_1}^{r_N} f(r) j_l(q r) r^2 dr
```
"""
@inline function bessel_transform(l::Int, r::AbstractVector{T},
                                  dr::Union{T,AbstractVector{T}}, f::AbstractVector{T},
                                  q::T)::T where {T<:Real}
    integrand = Vector{T}(undef, length(f))
    @inbounds @fastmath for i in eachindex(f)
        integrand[i] = f[i] * fast_sphericalbesselj(l, q * r[i])
    end
    return 4π * simpson(integrand, dr)
end

# Special case where l = 0 can use the (slightly) faster `fast_sphericalbesselj0`
@inline function bessel_transform(r::AbstractVector{T}, dr::Union{T,AbstractVector{T}},
                                  f::AbstractVector{T}, q::T)::T where {T<:Real}
    # This could be optimized for q = 0 where j₀(0) = 1.
    # It could also be optimized by substituting j₀ = sin(qr)/(qr) in the integrand
    # expression yielding `4π/q * integrator(r .* sin.(q .* r) .* f, dr)` where
    # the special cases of r = 0 and q = 0 would need to be properly handled.
    integrand = Vector{T}(undef, length(f))
    @inbounds @fastmath for i in eachindex(f)
        integrand[i] = f[i] * fast_sphericalbesselj0(q * r[i])
    end
    return 4π * simpson(integrand, dr)
end

@inline function bessel_transform(::Int, ::AbstractVector{T}, ::Union{T,AbstractVector{T}},
    ::Nothing, ::T)::Nothing where {T<:Real}
return nothing
end

@inline function bessel_transform(::AbstractVector{T}, ::Union{T,AbstractVector{T}},
                                  ::Nothing, ::T)::Nothing where {T<:Real}
    return nothing
end
