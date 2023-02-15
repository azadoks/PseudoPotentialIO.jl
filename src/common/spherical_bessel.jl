#TODO Implement the recursive algorithm to compare performance and accuracy
@doc raw"""
Spherical Bessel function of the first kind jₗ(x).
    
Consistent with https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
and with `SpecialFunctions.sphericalbesselj`.

Specialized for integer `l` in the range`0 <= l <= 5`.
"""
function fast_sphericalbesselj(l::Integer, x::T)::T where {T<:Real}
    l == 0 && return fast_sphericalbesselj0(x)
    l == 1 && return fast_sphericalbesselj1(x)
    l == 2 && return fast_sphericalbesselj2(x)
    l == 3 && return fast_sphericalbesselj3(x)
    l == 4 && return fast_sphericalbesselj4(x)
    l == 5 && return fast_sphericalbesselj5(x)
    return error("The case l = $l is not implemented")
end

function fast_sphericalbesselj0(x::T)::T where {T<:Real}
    iszero(x) && return one(T)
    return sin(x) / x
end

function fast_sphericalbesselj1(x::T)::T where {T<:Real}
    iszero(x) && return zero(T)
    sinx, cosx = sincos(x)
    invx2 = (1 / x)^2
    sinx_xcosx = fma(-cosx, x, sinx)
    return sinx_xcosx * invx2
end

function fast_sphericalbesselj2(x::T)::T where {T<:Real}
    iszero(x) && return zero(T)
    invx = 1 / x
    sinx, cosx = sincos(x)
    invx2 = invx * invx
    invxsinx = sinx * invx
    tmp1 = 3 * invx2
    tmp2 = tmp1 - 1
    return invxsinx * tmp2 - tmp1 * cosx
end

function fast_sphericalbesselj3(x::T)::T where {T<:Real}
    iszero(x) && return zero(T)
    invx = 1 / x
    sinx, cosx = sincos(x)
    invx2 = invx * invx
    xcosx = x * cosx
    invx4 = invx2 * invx2
    poly1 = -6 * invx2
    poly_cos = fma(-15, invx4, invx2)
    poly_sin = fma(15, invx4, poly1)
    return (sinx * poly_sin + xcosx * poly_cos)
end

function fast_sphericalbesselj4(x::T)::T where {T<:Real}
    iszero(x) && return zero(T)
    invx = 1 / x
    sinx, cosx = sincos(x)
    invx2 = invx * invx
    xsinx = x * sinx
    invx4 = invx2 * invx2
    invx6 = invx2 * invx4
    poly1 = fma(-45, invx4, invx2)
    poly2 = 10 * invx2
    poly_cos = fma(-105, invx4, poly2)
    poly_sin = fma(105, invx6, poly1)
    return (xsinx * poly_sin + cosx * poly_cos)
end

function fast_sphericalbesselj5(x::T)::T where {T<:Real}
    iszero(x) && return zero(T)
    invx = 1 / x
    sinx, cosx = sincos(x)
    invx2 = invx * invx
    xcosx = x * cosx
    invx4 = invx2 * invx2
    invx6 = invx2 * invx4
    poly1 = fma(105, invx4, -invx2)
    poly2 = 15 * invx2
    poly3 = fma(-420, invx4, poly2)
    poly_cos = fma(-945, invx6, poly1)
    poly_sin = fma(945, invx6, poly3)
    return (sinx * poly_sin + xcosx * poly_cos)
end