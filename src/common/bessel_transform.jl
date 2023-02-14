@doc raw"""
Bessel-Fourier (also Hankel) transform of order `l` of a function `f` on a radial mesh `r`.
The function `f` should be rapidly decaying to zero within the bounds of the mesh.

```math
4\pi \int_0^\inf f(r) j_l(q r) r^2 dr \approx
4\pi \int_{r_1}^{r_N} f(r) j_l(q r) r^2 dr
```
"""
function bessel_transform(l::Int, r::AbstractVector{T}, dr::Union{T,AbstractVector{T}},
                          f::AbstractVector{T}, q::T)::T where {T<:Real}
    integrand = _bessel_transform_integrand(l, r, f)
    return 4π * simpson(integrand, firstindex(f), lastindex(f), dr, q)
end

@inline function bessel_transform(::Int, ::AbstractVector{T}, ::Union{T,AbstractVector{T}},
                                  ::Nothing, ::T)::Nothing where {T<:Real}
    return nothing
end

@inline function bessel_transform(::AbstractVector{T}, ::Union{T,AbstractVector{T}},
                                  ::Nothing, ::T)::Nothing where {T<:Real}
    return nothing
end

@inbounds function _bessel_transform_integrand(l::Int, r::AbstractVector{T},
                                               f::AbstractVector{T})::Function where {T<:Real}
    fast_sphericalbesseljl = fast_sphericalbesselj_function(l)
    integrand(i::Int, q::T)::T = f[i] * fast_sphericalbesseljl(q * r[i])
    return integrand
end
