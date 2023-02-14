#TODO Implement the recursive algorithm to compare performance and accuracy
"""
Spherical Bessel function of the first kind jₗ(x).
    
Consistent with https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions
and with `SpecialFunctions.sphericalbesselj`.

Specialized for integer `l` in the range`0 <= l <= 5`.
"""
@fastmath function fast_sphericalbesselj(l::Integer, x::T)::T where {T}
    l == 0 && return fast_sphericalbesselj0(x)
    iszero(x) && return zero(T)
    l == 1 && return (sin(x) - cos(x) * x) / x^2
    l == 2 && return (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
    l == 3 && return (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
    l == 4 && return (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
    l == 5 &&
        return (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
    return error("The case l = $l is not implemented")
end

@inline @fastmath function fast_sphericalbesselj0(x::T)::T where {T}
    return iszero(x) ? one(T) : sin(x) / x
end

@inline @fastmath function fast_sphericalbesselj1(x::T)::T where {T}
    return iszero(x) ? x : (sin(x) - cos(x) * x) / x^2
end

@inline @fastmath function fast_sphericalbesselj2(x::T)::T where {T}
    return iszero(x) ? x : (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
end

@inline @fastmath function fast_sphericalbesselj3(x::T)::T where {T}
    return iszero(x) ? x : (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
end

@inline @fastmath function fast_sphericalbesselj4(x::T)::T where {T}
    return iszero(x) ? x : (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
end

@inline @fastmath function fast_sphericalbesselj5(x::T)::T where {T}
    return iszero(x) ? x :
           (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
end

"""
Get the spherical Bessel function of the first kind jₗ(x) for a given `l`.
"""
function fast_sphericalbesselj_function(l::Integer)::Function
    l == 0 && return fast_sphericalbesselj0
    l == 1 && return fast_sphericalbesselj1
    l == 2 && return fast_sphericalbesselj2
    l == 3 && return fast_sphericalbesselj3
    l == 4 && return fast_sphericalbesselj4
    l == 5 && return fast_sphericalbesselj5
    return error("The case l = $l is not implemented")
end