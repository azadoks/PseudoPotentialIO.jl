@doc raw"""
Trapezoidal rule integration for a function `f(x)` on a grid with grid spacing `dx`.

For a non-uniform grid:
```math
\int_a^b f(x) dx \approx \sum_{k=2}^N \frac{f(x_{k-1}) + f(x_k)}{2} (x_k - x_{k-1})
```

For a uniform grid:
```math
\int_a^b f(x) dx \approx
\Delta x \left( \sum_{k=2}^{N-1} f(x_k) + \frac{f(x_1) + f(x_N)}{2} \right)
```
"""
@inline function trapezoid(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    s = zero(T)
    @inbounds @simd for i in (firstindex(f) + 1):lastindex(f)
        s += (f[i] + f[i - 1]) * dx[i]
    end
    return s / 2
end

@inline function trapezoid(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = (first(f) + last(f)) / 2
    s += sum(@view f[begin + 1:end - 1])
    return s * dx
end

@doc raw"""
Simpson's rule integration for a function `f(x)` on a grid with grid spacing `dx`.

Performs better than the trapezoidal rule on logarithmic grids.

For a uniform grid with an odd number of grid points:
```math
\int_a^b f(x) dx \approx
\frac{\Delta x}{3} \left[
f(x_1) +
4 \left( \sum_{i=2, i \mathrm{even}}^{N-1} f(x_i) \right) +
2 \left( \sum_{i=3, i \mathrm{odd}}^{N-2} f(x_i) +
f(x_N) \right) \right]
```
"""
@inline function simpson(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    s = f[begin] * dx[begin]
    s += @inbounds sum(i -> 4 * f[i] * dx[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += @inbounds sum(i -> 2 * f[i] * dx[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    @inbounds s += length(f) % 2 == 1 ? f[end] * dx[end] : -f[end - 1] * dx[end - 1]
    return s / 3
end

@inline function simpson(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = f[begin]
    s += @inbounds 4 * sum(i -> f[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += @inbounds 2 * sum(i -> f[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    @inbounds s += length(f) % 2 == 1 ? f[end] : -f[end - 1]
    return s / 3 * dx
end
