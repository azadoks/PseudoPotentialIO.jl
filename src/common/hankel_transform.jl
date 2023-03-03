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
function hankel_transform(::OrbitalLike, l::Int, r::AbstractVector,
                          dr::Union{Real,AbstractVector}, f::AbstractVector;
                          i_start=firstindex(f), i_stop=lastindex(f))
    jₗ = fast_sphericalbesselj(l)
    f̃(q) = 4π * simpson(i -> f[i] * jₗ(q * r[i]), i_start, i_stop, dr)
    f̃(Q::AbstractVector) = f̃(norm(Q))
    return f̃
end

function hankel_transform(::DensityLike, r::AbstractVector, dr::Union{Real,AbstractVector},
                          f::AbstractVector; i_start=firstindex(f), i_stop=lastindex(f))
    function f̃(q)
        integrand(i) = f[i] * fast_sphericalbesselj0(q * r[i])
        return 4π * simpson(integrand, i_start, i_stop, dr)
    end
    f̃(Q::AbstractVector) = f̃(norm(Q))
    return f̃
end

@inline function hankel_transform(::OrbitalLike, ::Int, ::AbstractVector,
                                  ::Union{Real,AbstractVector}, ::Nothing; i_start=nothing,
                                  i_stop=nothing)
    return _ -> nothing
end

@inline function hankel_transform(::DensityLike, ::AbstractVector,
                                  ::Union{Real,AbstractVector}, ::Nothing; i_start=nothing,
                                  i_stop=nothing)
    return _ -> nothing
end
