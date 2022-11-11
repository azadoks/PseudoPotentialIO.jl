"""
$(SIGNATURES)

Trapezoidal rule integration for a function `f(x)` on a grid with grid spacing `dx`.
"""
@inline function trapezoid(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    s = zero(T)
    @inbounds @simd for i in eachindex(f)
        s += f[i] * dx[i]
    end
    return s
end

@inline function trapezoid(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    return sum(f) * dx
end

"""
$(SIGNATURES)

Simpson's rule integration for a function `f(x)` on a grid with grid spacing `dx`.
"""
@inline function simpson(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    n = length(f)
    s = f[begin] * dx[begin]
    s += @inbounds sum(i -> 4 * f[i] * dx[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += @inbounds sum(i -> 2 * f[i] * dx[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    @inbounds s += n % 2 == 1 ? f[n] * dx[n] : -f[n - 1] * dx[n - 1]
    return s / T(3)
end

@inline function simpson(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = f[begin]
    s += @inbounds sum(i -> 4 * f[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += @inbounds sum(i -> 2 * f[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    @inbounds s += length(f) % 2 == 1 ? f[end] : -f[end - 1]
    return s / T(3) * dx
end
