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
trapezoid
@inbounds function trapezoid(f::AbstractVector{T},
                             dx::AbstractVector{T})::T where {T<:Real}
    s = zero(T)
    @simd for i in (firstindex(f) + 1):lastindex(f)
        s += (f[i] + f[i - 1]) * dx[i]
    end
    return s / 2
end

@inbounds function trapezoid(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = (first(f) + last(f)) / 2
    s += sum(@view f[(begin + 1):(end - 1)])
    return s * dx
end

function trapezoid(f::Function, i_start::Int, i_stop::Int, dx::AbstractVector{T},
                   q::T)::T where {T<:Real}
    s = zero(T)
    fa = f(i_start, q)
    for i in (i_start + 1):i_stop
        fb = f(i, q)
        s += (fb + fa) * dx[i]
        fa = fb
    end
    return s / 2
end

function trapezoid(f::Function, i_start::Int, i_stop::Int, dx::T, q::T)::T where {T<:Real}
    s = (f(i_start, q) + f(i_stop)) / 2
    for i in (i_start + 1):i_stop
        s += f(i, q)
    end
    return s * dx
end

@inbounds function simpson(f::AbstractVector{T}, dx::AbstractVector{T})::T where {T<:Real}
    s = f[begin] * dx[begin]
    s += sum(i -> 4 * f[i] * dx[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += sum(i -> 2 * f[i] * dx[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    s += length(f) % 2 == 1 ? f[end] * dx[end] : -f[end - 1] * dx[end - 1]
    return s / 3
end

@inbounds function simpson(f::AbstractVector{T}, dx::T)::T where {T<:Real}
    s = f[begin]
    s += 4 * sum(i -> f[i], (firstindex(f) + 1):2:(lastindex(f) - 1))
    s += 2 * sum(i -> f[i], (firstindex(f) + 2):2:(lastindex(f) - 1))
    s += length(f) % 2 == 1 ? f[end] : -f[end - 1]
    return s / 3 * dx
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
simpson
@inbounds function simpson(f::Function, i_start::Int, i_stop::Int, dx::AbstractVector{T},
                           q::T)::T where {T<:Real}
    s = f(i_start, q) * dx[i_start]
    s += sum(i -> 4 * f(i, q) * dx[i], (i_start + 1):2:(i_stop - 1))
    s += sum(i -> 2 * f(i, q) * dx[i], (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start) % 2 == 1 ? f(i_stop, q) * dx[i_stop] :
         -f(i_stop - 1, q) * dx[i_stop - 1]
    return s / 3
end

@inbounds function simpson(f::Function, i_start::Int, i_stop::Int, dx::T,
                           q::T)::T where {T<:Real}
    s = f(i_start, q)
    s += 4 * sum(i -> f(i, q), (i_start + 1):2:(i_stop - 1))
    s += 2 * sum(i -> f(i, q), (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start) % 2 == 1 ? f(i_stop, q) : -f(i_stop - 1, q)
    return s / 3 * dx
end
