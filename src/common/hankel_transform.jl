@doc raw"""
Hankel / Bessel-Fourier transform of order `l` of a function `f` on a radial mesh `r`.
The function `f` should be rapidly decaying to zero within the bounds of the mesh.

The radial transform is defined as:

```math
4\pi \int_0^{\infty} f(r) j_l(q r) r^2 dr \approx
4\pi \int_{r_1}^{r_N} f(r) j_l(q r) r^2 dr
```

where ``j_l(x)`` is the spherical Bessel function of the first kind at order ``l``.
"""
hankel_transform
@inbounds function hankel_transform(f::AbstractVector, l::Int, r::AbstractVector,
                                    dr::Union{Real,AbstractVector}; i_start=firstindex(f),
                                    i_stop=lastindex(f), integrator=simpson)
    jₗ = fast_sphericalbesselj(l)
    @inbounds f̃(q) = 4π * integrator(i -> f[i] * jₗ(q * r[i]), i_start, i_stop, dr)
    @inbounds f̃(Q::AbstractVector) = f̃(norm(Q))
    return f̃
end

@inline function hankel_transform(::Nothing, ::Int, ::AbstractVector,
                                  ::Union{Real,AbstractVector}; i_start=nothing,
                                  i_stop=nothing)
    return _ -> nothing
end
