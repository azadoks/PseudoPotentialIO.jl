abstract type BesselTransformQuantityType end
struct OrbitalLike <: BesselTransformQuantityType end
struct DensityLike <: BesselTransformQuantityType end

@doc raw"""
Bessel-Fourier (also Hankel) transform of order `l` of a function `f` on a radial mesh `r`.
The function `f` should be rapidly decaying to zero within the bounds of the mesh.

```math
4\pi \int_0^\inf f(r) j_l(q r) r^2 dr \approx
4\pi \int_{r_1}^{r_N} f(r) j_l(q r) r^2 dr
```
"""
function bessel_transform(quantity_type::BesselTransformQuantityType, l::Int,
                          r::AbstractVector, dr::Union{Real,AbstractVector},
                          f::AbstractVector, q)
    integrand = _bessel_transform_integrand(quantity_type, l, r, f, q)
    return 4π * simpson(integrand, firstindex(f), lastindex(f), dr)
end

@inline function bessel_transform(::BesselTransformQuantityType, ::Int, ::AbstractVector,
                                  ::Union{Real,AbstractVector}, ::Nothing,
                                  _)::Nothing
    return nothing
end

@inline function bessel_transform(::BesselTransformQuantityType, ::AbstractVector,
                                  ::Union{Real,AbstractVector},
                                  ::Nothing, _)::Nothing
    return nothing
end

@inbounds function _bessel_transform_integrand(::OrbitalLike, l::Int, r::AbstractVector,
                                               f::AbstractVector, q)::Function
    integrand(i::Int) = f[i] * fast_sphericalbesselj(l, q * r[i])
    return integrand
end

@inbounds function _bessel_transform_integrand(::DensityLike, l::Int, r::AbstractVector,
                                               f::AbstractVector, q)::Function
    integrand(i::Int) = r[i]^2 * f[i] * fast_sphericalbesselj(l, q * r[i])
    return integrand
end