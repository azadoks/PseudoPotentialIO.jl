@doc raw"""
Hankel / Bessel-Fourier transform of order `l` of a function `f` on a radial mesh `r`.
The function `f` should be rapidly decaying to zero within the bounds of the mesh.
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

function hankel_transform(r::AbstractVector, r²f::AbstractVector,
                          dr::AbstractVector, Δr::Union{Real,AbstractVector}, l::Int;
                          quadrature_method=Trapezoid())
    jₗ = fast_sphericalbesselj(l)
    @inbounds function f̃(q::T) where {T<:Real}
        return 4T(π) *
               integrate(eachindex(r²f), i -> r²f[i] * jₗ(q * r[i]), dr, Δr, quadrature_method)
    end
    @inbounds f̃(Q::AbstractVector) = f̃(norm(Q))
    return f̃
end

function hankel_transform(r::AbstractVector, r²f, dr::AbstractVector, Δr::Union{Real,AbstractVector}, l::Int;
                          quadrature_method=Trapezoid())
    jₗ = fast_sphericalbesselj(l)
    @inbounds function f̃(q::T) where {T<:Real}
        return 4T(π) * integrate(r, ri -> r²f(ri) * jₗ(q * ri), dr, Δr, quadrature_method)
    end
    @inbounds f̃(Q::AbstractVector) = f̃(norm(Q))
    return f̃
end

function hankel_transform(::AbstractVector, ::Nothing, ::Union{Real,AbstractVector}, ::Int;
                          quadrature_method=Trapezoid())
    return _ -> nothing
end

function hankel_transform!(f̃::AbstractVector, q::AbstractVector{T}, r::AbstractVector,
                           r²f::AbstractVector, dr::AbstractVector, Δr::Union{Real,AbstractVector}, l::Int;
                           quadrature_method=Trapezoid()) where {T<:Real}
    jₗ = fast_sphericalbesselj(l)
    map!(f̃, q) do qj
        return 4T(π) *
               integrate(eachindex(r²f), i -> r²f[i] * jₗ(qj * r[i]), dr, Δr, quadrature_method)
    end
    return f̃
end

function hankel_transform!(f̃::AbstractVector, q::AbstractVector{T}, r::AbstractVector, r²f,
                           dr::AbstractVector, Δr::Union{Real,AbstractVector}, l::Int;
                           quadrature_method=Trapezoid()) where {T<:Real}
    jₗ = fast_sphericalbesselj(l)
    map!(f̃, q) do qj
        return 4T(π) * integrate(r, ri -> r²f(ri) * jₗ(qj * ri), dr, Δr, quadrature_method)
    end
    return f̃
end

function hankel_transform!(f̃::AbstractVector, ::AbstractVector, ::AbstractVector,
                           ::Nothing, ::Union{Real,AbstractVector}, ::Int;
                           quadrature_method=Trapezoid())
    f̃ .= nothing
    return f̃
end
