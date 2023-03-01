abstract type HankelTransformQuantityType end
struct OrbitalLike <: HankelTransformQuantityType end
struct DensityLike <: HankelTransformQuantityType end

@doc raw"""
Hankel / Bessel-Fourier transform of order `l` of a function `f` on a radial mesh `r`.
The function `f` should be rapidly decaying to zero within the bounds of the mesh.

```math
4\pi \int_0^\inf f(r) j_l(q r) r^2 dr \approx
4\pi \int_{r_1}^{r_N} f(r) j_l(q r) r^2 dr
```
"""
function hankel_transform(quantity_type::HankelTransformQuantityType, l::Int,
                          r::AbstractVector, dr::Union{Real,AbstractVector},
                          f::AbstractVector, q)
    integrand = _hankel_transform_integrand(quantity_type, l, r, f, q)
    return 4π * simpson(integrand, firstindex(f), lastindex(f), dr)
end

@inline function hankel_transform(::HankelTransformQuantityType, ::Int, ::AbstractVector,
                                  ::Union{Real,AbstractVector}, ::Nothing,
                                  _)::Nothing
    return nothing
end

@inline function hankel_transform(::HankelTransformQuantityType, ::AbstractVector,
                                  ::Union{Real,AbstractVector},
                                  ::Nothing, _)::Nothing
    return nothing
end

@inbounds function _hankel_transform_integrand(::OrbitalLike, l::Int, r::AbstractVector,
                                               f::AbstractVector, q)::Function
    integrand(i::Int) = f[i] * fast_sphericalbesselj(l, q * r[i])
    return integrand
end

@inbounds function _hankel_transform_integrand(::DensityLike, l::Int, r::AbstractVector,
                                               f::AbstractVector, q)::Function
    integrand(i::Int) = r[i]^2 * f[i] * fast_sphericalbesselj(l, q * r[i])
    return integrand
end

function hankel_transform_function(::OrbitalLike, l::Int,
                                   r::AbstractVector, dr::Union{Real,AbstractVector},
                                   f::AbstractVector)
    j = fast_sphericalbesselj_function(l)
    hankel(q) = 4π * simpson(i -> f[i] * j(q * r[i]), firstindex(f), lastindex(f), dr)
    return hankel
end

function hankel_transform_function(::DensityLike, r::AbstractVector,
                                   dr::Union{Real,AbstractVector}, f::AbstractVector)
    function hankel(q)
        integrand(i) = r[i]^2 * f[i] * fast_sphericalbesselj0(q * r[i])
        return 4π * simpson(integrand, firstindex(f), lastindex(f), dr)
    end
    return hankel
end

@inline function hankel_transform_function(::HankelTransformQuantityType, ::Int,
                                           ::AbstractVector, ::Union{Real,AbstractVector},
                                           ::Nothing)
    return _ -> nothing
end

@inline function hankel_transform_function(::HankelTransformQuantityType, ::AbstractVector,
                                           ::Union{Real,AbstractVector}, ::Nothing)
    return _ -> nothing
end
