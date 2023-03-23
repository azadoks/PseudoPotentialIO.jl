@doc raw"""
Simpson's rule integration for a function `f(x)` on a grid with grid spacing `dx`.
Performs better than the trapezoidal rule on logarithmic grids.

For a uniform grid with an odd number of grid points:
```math
\int_a^b f(x) dx \approx
\frac{\Delta x}{3} \left[
f(x_1) +
4 \left( \sum_{i=2, i_\mathrm{even}}^{N-1} f(x_i) \right) +
2 \left( \sum_{i=3, i_\mathrm{odd}}^{N-2} f(x_i) +
f(x_N) \right) \right]
```
"""
simpson
@inbounds function simpson(f, i_start::Int, i_stop::Int, dx::AbstractVector)
    s = f(i_start) * dx[i_start]
    s += sum(i -> 4 * f(i) * dx[i], (i_start + 1):2:(i_stop - 1))
    s += sum(i -> 2 * f(i) * dx[i], (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start + 1) % 2 == 1 ? f(i_stop) * dx[i_stop] :
         -f(i_stop - 1) * dx[i_stop - 1]
    return s / 3
end

@inbounds function simpson(f, i_start::Int, i_stop::Int, dx)
    s = f(i_start)
    s += 4 * sum(i -> f(i), (i_start + 1):2:(i_stop - 1))
    s += 2 * sum(i -> f(i), (i_start + 2):2:(i_stop - 1))
    s += (i_stop - i_start + 1) % 2 == 1 ? f(i_stop) : -f(i_stop - 1)
    return s / 3 * dx
end

@inbounds function rectangle(f, i_start::Int, i_stop::Int, dx::AbstractVector)
    return sum(i -> f(i) * dx[i], i_start:i_stop)
end

@inbounds function rectangle(f, i_start::Int, i_stop::Int, dx)
    return sum(f, i_start:i_stop) * dx
end

@inbounds function trapezoid(f, i_start::Int, i_stop::Int, dx::AbstractVector)
    return sum(i -> (f(i) + f(i + 1)) * dx[i], i_start:(i_stop - 1)) / 2
end

@inbounds function trapezoid(f, i_start::Int, i_stop::Int, dx)
    return dx * (sum(f, (i_start + 1):(i_stop - 1)) + (f(i_start) + f(i_stop)) / 2)
end

# Same as ABINIT `ctrap`
# https://github.com/abinit/abinit/blob/6e17fde483be5e66beaa09debf3f0f6dcf4d98cb/shared/common/src/28_numeric_noabirule/m_numeric_tools.F90#L2728
function abinit_corrected_trapezoid(f, i_start::Int, i_stop::Int, dx::Real)
    n = i_stop - i_start + 1
    if n >= 10
        endpoint = (
            23.75(f(i_start    ) + f(i_stop    )) +
            95.10(f(i_start + 1) + f(i_stop - 1)) +
            55.20(f(i_start + 2) + f(i_stop - 2)) +
            79.30(f(i_start + 3) + f(i_stop - 3)) +
            70.65(f(i_start + 4) + f(i_stop - 4))
        ) / 72
        return (sum(f, i_start+6:i_stop-5) + endpoint) * dx
    elseif n >= 8
        endpoint = (
            17(f(i_start    ) + f(i_stop    )) +
            59(f(i_start + 1) + f(i_stop - 1)) +
            43(f(i_start + 2) + f(i_stop - 2)) +
            49(f(i_start + 3) + f(i_stop - 3))
        ) / 48
        return n == 8 ? endpoint * dx : (endpoint + f(i_start + 4)) * dx
    elseif n == 7
        return (
            17(f(i_start    ) + f(i_stop)    ) +
            59(f(i_start + 1) + f(i_stop - 1)) +
            43(f(i_start + 2) + f(i_stop - 2)) +
            50(f(i_start + 3)                )
        ) / 48 * dx
    elseif n == 6
        return (
            17(f(i_start    ) + f(i_stop    )) +
            59(f(i_start + 1) + f(i_stop - 1)) +
            44(f(i_start + 2) + f(i_stop - 2))
        ) / 48 * dx
    elseif n == 5
        return (
             (f(i_start    ) + f(i_start + 4)) +
            4(f(i_start + 1) + f(i_start + 3)) +
            2(f(i_start + 2                 ))
        ) / 3 * dx
    elseif n == 4
        return (
            3(f(i_start   ) + f(i_start + 1)) +
            9(f(i_stop - 1) + f(i_stop     ))
        ) / 8 * dx
    elseif n == 3
        return (
             (f(i_start    ) + f(i_start + 2)) +
            4(f(i_start + 1) + f(i_start + 1))
        ) / 3 * dx
    elseif n == 2
        return (f(i_start) + f(i_stop)) / 2 * dx
    elseif n == 1
        return f(i_start) * dx
    end
end

# Same as q-e `simpson`
# https://github.com/QEF/q-e/blob/48b24e82928af44ba63d04f37c99b0224cbc506d/upflib/simpsn.f90#L9
function qe_simpson(f, i_start::Int, i_stop::Int, dx::AbstractVector)
    n = i_stop - i_start + 1
    s = sum((i_start + 1):(i_stop - 1)) do i
        2abs(mod(i, 2)-2) * f(i) * dx[i]
    end
    if mod(n, 2) == 1
        s = (s + f(i_start) * dx[i_start] + f(i_stop) * dx[i_stop]) / 3
    else
        s = (s + f(i_start) * dx[i_start] + f(i_stop-1) * dx[i_stop-1]) / 3
    end
    return s
end

function qe_simpson(f, i_start::Int, i_stop::Int, dx)
    n = i_stop - i_start + 1
    s = sum((i_start + 1):(i_stop - 1)) do i
        2abs(mod(i, 2)-2) * f(i) * dx
    end
    if mod(n, 2) == 1
        s = (s + f(i_start) * dx + f(i_stop) * dx) / 3
    else
        s = (s + f(i_start) * dx + f(i_stop-1) * dx) / 3
    end
    return s
end

# Same as q-e `simpson_cp90`
# https://github.com/QEF/q-e/blob/48b24e82928af44ba63d04f37c99b0224cbc506d/upflib/simpsn.f90#L61
function cp90_simpson(f, i_start::Int, i_stop::Int, dx::AbstractVector)
    n = i_stop - i_start + 1
    n < 8 && throw(ArgumentError("Minimum 8 points are required"))

    s = (
        (f(i_start    ) * dx[i_start    ] + f(i_stop    ) * dx[i_stop    ]) *  109 +
        (f(i_start + 1) * dx[i_start + 1] + f(i_stop - 1) * dx[i_stop - 1]) * -5   +
        (f(i_start + 2) * dx[i_start + 2] + f(i_stop - 2) * dx[i_stop - 2]) *  63  +
        (f(i_start + 3) * dx[i_start + 3] + f(i_stop - 3) * dx[i_stop - 3]) *  49 
    ) / 48 
    for i in (i_start + 4):(i_stop - 4)
        s += f(i) * dx[i]
    end
    return s
end

function cp90_simpson(f, i_start::Int, i_stop::Int, dx)
    n = i_stop - i_start + 1
    n < 8 && throw(ArgumentError("Minimum 8 points are required"))

    s = (
        (f(i_start    ) + f(i_stop    )) *  109 +
        (f(i_start + 1) + f(i_stop - 1)) * -5   +
        (f(i_start + 2) + f(i_stop - 2)) *  63  +
        (f(i_start + 3) + f(i_stop - 3)) *  49 
    ) * dx / 48
    for i in (i_start + 4):(i_stop - 4)
        s += f(i) * dx
    end
    return s
end
