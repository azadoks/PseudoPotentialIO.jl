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
@inbounds function simpson(f::Function, i_start::Int, i_stop::Int,
                           dx::AbstractVector{T})::T where {T<:Real}
    s = f(i_start) * dx[i_start]
    s += sum(i -> 4 * f(i) * dx[i], (i_start + 1):2:(i_stop - 1))
    s += sum(i -> 2 * f(i) * dx[i], (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start + 1) % 2 == 1 ? f(i_stop) * dx[i_stop] :
         -f(i_stop - 1) * dx[i_stop - 1]
    return s / 3
end

@inbounds function simpson(f::Function, i_start::Int, i_stop::Int, dx::T)::T where {T<:Real}
    s = f(i_start)
    s += 4 * sum(i -> f(i), (i_start + 1):2:(i_stop - 1))
    s += 2 * sum(i -> f(i), (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start + 1) % 2 == 1 ? f(i_stop) : -f(i_stop - 1)
    return s / 3 * dx
end
