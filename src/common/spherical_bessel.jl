#TODO Implement the recursive algorithm to compare performance and accuracy
#TODO Try to remove branching in `fast_sphericalbesselj`
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

function fast_sphericalbesselj_function(l::Integer)
    l == 0 && return fast_sphericalbesselj0
    l == 1 && return fast_sphericalbesselj1
    l == 2 && return fast_sphericalbesselj2
    l == 3 && return fast_sphericalbesselj3
    l == 4 && return fast_sphericalbesselj4
    l == 5 && return fast_sphericalbesselj5
    return error("The case l = $l is not implemented")
end

function fast_sphericalbesselj0(x::T)::T where {T<:Real}
    if !iszero(x)
        return sin(x) / x
    end
    return one(T)
end

function fast_sphericalbesselj1(x::T)::T where {T<:Real}
    if !iszero(x)
        # (sin(x) - cos(x) * x) / x^2
        sinx, cosx = sincos(x)
        invx2 = (1 / x)^2
        sinx_xcosx = fma(-cosx, x, sinx)
        return sinx_xcosx * invx2
    end
    return zero(T)
end

function fast_sphericalbesselj2(x::T)::T where {T<:Real}
    if !iszero(x)
        # (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
        invx = 1 / x
        sinx, cosx = sincos(x)
        invx2 = invx * invx
        invxsinx = sinx * invx
        tmp1 = 3 * invx2
        tmp2 = tmp1 - 1
        return invxsinx * tmp2 - tmp1 * cosx
    end
    return zero(T)
end

function fast_sphericalbesselj3(x::T)::T where {T<:Real}
    if !iszero(x)
        # (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
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
    return zero(T)
end

function fast_sphericalbesselj4(x::T)::T where {T<:Real}
    if !iszero(x)
        # (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
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
    return zero(T)
end

function fast_sphericalbesselj5(x::T)::T where {T<:Real}
    if !iszero(x)
        # (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
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
    return zero(T)
end