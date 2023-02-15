#TODO Implement the recursive algorithm to compare performance and accuracy
@doc raw"""
Spherical Bessel function of the first kind jₗ(x).
    
Consistent with https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
and with `SpecialFunctions.sphericalbesselj`.

Specialized for integer `l` in the range`0 <= l <= 5`.
"""
function fast_sphericalbesselj(l::Integer, x::T)::T where {T<:Real}
    l == 0 && return fast_sphericalbesselj0(x)
    iszero(x) && return zero(T)
    l == 1 && return _sphericalbesselj1(x)
    # (sin(x) - cos(x) * x) / x^2
    l == 2 && return _sphericalbesselj2(x)
    # (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
    l == 3 && return _sphericalbesselj3(x)
    # (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
    l == 4 && return _sphericalbesselj4(x)
    # (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
    l == 5 && return _sphericalbesselj5(x)
    # return (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
    return error("The case l = $l is not implemented")
end

@inline function fast_sphericalbesselj0(x::T)::T where {T<:Real}
    return iszero(x) ? one(T) : sin(x) / x
end

@inline function _sphericalbesselj1(x::T)::T where {T<:Real}
    sinx = sin(x)
    cosx = cos(x)
    invx2 = (1 / x)^2
    sinx_xcosx = fma(-cosx, x, sinx)
    return sinx_xcosx * invx2
end

@inline function _sphericalbesselj2(x::T)::T where {T<:Real}
    invx = 1 / x
    sinx = sin(x)
    cosx = cos(x)
    invx2 = invx^2
    xsinx = x * sinx
    invx4 = invx2^2
    poly_cos = -3 * invx2
    poly_sin = fma(3, invx4, -invx2)
    return xsinx * poly_sin + cosx * poly_cos
end

@inline function _sphericalbesselj3(x::T)::T where {T<:Real}
    invx = 1 / x
    sinx = sin(x)
    cosx = cos(x)
    invx2 = invx^2
    xcosx = x * cosx
    invx4 = invx2^2
    poly1 = -6 * invx2
    poly_cos = fma(-15, invx4, invx2)
    poly_sin = fma(15, invx4, poly1)
    return (sinx * poly_sin + xcosx * poly_cos)
end

@inline function _sphericalbesselj4(x::T)::T where {T<:Real}
    invx = 1 / x
    sinx = sin(x)
    cosx = cos(x)
    invx2 = invx^2
    xsinx = x * sinx
    invx4 = invx2^2
    invx6 = invx2 * invx4
    poly1 = fma(-45, invx4, invx2)
    poly2 = 10 * invx2
    poly_cos = fma(-105, invx4, poly2)
    poly_sin = fma(105, invx6, poly1)
    return (xsinx * poly_sin + cosx * poly_cos)
end

@inline function _sphericalbesselj5(x::T)::T where {T<:Real}
    invx = 1 / x
    sinx = sin(x)
    cosx = cos(x)
    invx2 = invx^2
    xcosx = x * cosx
    invx4 = invx2^2
    invx6 = invx2 * invx4
    poly1 = fma(105, invx4, -invx2)
    poly2 = 15 * invx2
    poly3 = fma(-420, invx4, poly2)
    poly_cos = fma(-945, invx6, poly1)
    poly_sin = fma(945, invx6, poly3)
    return (sinx * poly_sin + xcosx * poly_cos)
end