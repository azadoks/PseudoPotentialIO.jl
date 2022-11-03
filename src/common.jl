"""
$(SIGNATURES)

Compute the value of a linear mesh `a * i + b` at index `i`.
"""
linear_mesh(i::Int, a::T, b::T) where {T<:Real} = a * i + b

"""
$(SIGNATURES)

Compute the value of a logarithmic mesh `b * exp(a * (i - 1))` at index `i`.
"""
logarithmic_mesh1(i::Int, a::T, b::T) where {T<:Real} = b * exp(a * (i - 1))
function logarithmic_mesh1(i::Int, xmin::T, dx::T, z::T) where {T<:Real}
    return exp(xmin) * exp((i - 1) * dx) / z
end

"""
$(SIGNATURES)

Compute the value of a logarithmic mesh `b * (exp(a * (i - 1)) - 1)` at index `i`.
"""
logarithmic_mesh2(i::Int, a::T, b::T) where {T<:Real} = b * (exp((i - 1) * a) - 1)
function logarithmic_mesh2(i::Int, xmin::T, dx::T, z::T) where {T<:Real}
    return exp(xmin) * (exp((i - 1) * dx) - 1) / z
end

"""
$(SIGNATURES)

Guess whether a numerical mesh is linear or one of two types of logarithmic mesh used in 
UPF pseudopotentials.
"""
function guess_mesh_type(r::Vector{T}, rab::Vector{T}) where {T<:Real}
    nr = length(r)
    # Try linear
    a = r[2] - r[1]
    b = r[1] - a
    rguess = linear_mesh.(1:nr, a, b)
    if all(rguess .≈ r) & all(round.(rab, digits=4) .≈ a)
        return ("linear", a, b)
    end
    # Try log1
    a = log(r[2] / r[1])  # dx
    b = r[2] / exp(a)  # exp(xmin) / zmesh
    rguess = logarithmic_mesh1.(1:nr, a, b)
    if all(rguess .≈ r) && all(rab .≈ a .* r)
        return ("log_1", a, b)
    end
    # Try log2
    b = (r[2]^2 - r[3] * r[1]) / (r[1] + r[3] - 2 * r[2])
    a = log((r[2] + b) / (r[1] + b))
    rguess = logarithmic_mesh2.(1:nr, a, b)
    if all(rguess .≈ r) && all(rab .≈ a .* r .+ a * b)
        return ("log_2", a, b)
    end
    return ("unknown", NaN, NaN)
end

"""
$(SIGNATURES)

Returns the spherical Bessel function of the first kind j_l(x). Consistent with
https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions and with
`SpecialFunctions.sphericalbesselj`. Specialized for integer `0 <= l <= 5`.
"""
@fastmath function sphericalbesselj_fast(l::Integer, x::T)::T where {T}
    if l == 0
        iszero(x) && return one(T)
        return sin(x) / x
    end

    iszero(x) && return zero(T)

    l == 1 && return (sin(x) - cos(x) * x) / x^2
    l == 2 && return (sin(x) * (3 - x^2) + cos(x) * (-3x)) / x^3
    l == 3 && return (sin(x) * (15 - 6x^2) + cos(x) * (x^3 - 15x)) / x^4
    l == 4 && return (sin(x) * (105 - 45x^2 + x^4) + cos(x) * (10x^3 - 105x)) / x^5
    l == 5 &&
        return (sin(x) * (945 - 420x^2 + 15x^4) + cos(x) * (-945x + 105x^3 - x^5)) / x^6
    return error("The case l = $l is not implemented")
end

"""
$(SIGNATURES)

Trapezoidal rule integration for a function `f(x)` on a grid with grid spacing `dx`.
"""
@inline function trapezoid(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    return dot(f, dx)
    # s = zero(T)
    # @inbounds for i in eachindex(f)
    # 	s += f[i] * dx[i]
    # end
    # return s
end

@inline function trapezoid(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    return sum(f) * dx
    # s = zero(T)
    # @inbounds for i in eachindex(f)
    # 	s += f[i] * dx[i]
    # end
    # return s
end

"""
$(SIGNATURES)

Simpson's rule integration for a function `f(x)` on a grid with grid spacing `dx`.
"""
@inline function simpson(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    s = f[begin] * dx[begin]
    s += sum(i -> 4 * f[i] * dx[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += sum(i -> 2 * f[i] * dx[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    s += length(f) % 2 == 1 ? f[end] * dx[end] : -f[end - 1] * dx[end - 1]
    return s / T(3)
    # for i in (firstindex(f) + 1):2:(lastindex(f) - 1)
    # 	s += 4 * f[i] * dx[i]
    # end
    # for i in (firstindex(f) + 2):2:(lastindex(f - 1))
    # 	s += 2 * f[i] * dx[i]
    # end 
end

@inline function simpson(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = f[begin]
    s += sum(i -> 4 * f[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += sum(i -> 2 * f[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    s += length(f) % 2 == 1 ? f[end] : -f[end - 1]
    return s / T(3) * dx
end
