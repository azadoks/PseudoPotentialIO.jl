function integrate(::AbstractVector, y::AbstractVector, ::AbstractVector, Δx::Real, ::Trapezoid)
    length(y) >= 2 || error("x and y must contain at least two elements")

    return (sum(y) - (y[begin] + y[end]) / 2) * Δx
end

function integrate(x::AbstractVector, y::AbstractVector, ::AbstractVector, Δx::AbstractVector, ::Trapezoid)
    length(y) >= 2 || error("x and y must contain at least two elements")

    ifirst = firstindex(x)
    ilast = lastindex(x)

    integral = (y[begin] + y[begin + 1]) * Δx[begin]
    for i in (ifirst + 1):(ilast - 1)
        integral += (y[i] + y[i + 1]) * Δx[i]
    end
    return integral / 2
end

function integrate(x::AbstractVector, y::AbstractVector, ::AbstractVector, Δx::Real, ::Simpson)
    length(y) >= 4 || error("x and y must contain at least four elements")

    N = length(x) - 1  # Number of intervals
    ifirst = firstindex(y)
    ilast = lastindex(y)

    if !isodd(N)
        integral = y[begin]
        integral += 4 * sum(y[(ifirst + 1):2:(ilast - 1)])
        integral += 2 * sum(y[(ifirst + 2):2:(ilast - 1)])
        integral += y[end]
        return integral / 3 * Δx
    else
        integral = 5 / 4 * y[begin]
        integral += 13 / 4 * y[begin + 1]
        integral += 3 * sum(y[(ifirst + 2):(ilast - 2)])
        integral += 13 / 4 * y[end - 1]
        integral += 5 / 4 * y[end - 1]
        return integral / 3 * Δx
    end
end

function integrate(x::AbstractVector, y::AbstractVector, ::AbstractVector, Δx::AbstractVector, ::Simpson)
    length(y) >= 4 || error("y must contain at least four elements")

    N = length(x) - 1  # Number of intervals
    istop = isodd(N) ? lastindex(x) - 3 : lastindex(x) - 2

    integral = 0
    for i in firstindex(x):2:istop
        prefac = (Δx[i] + Δx[i + 1]) / 6
        integral += y[i] * prefac * (2 - Δx[i + 1] / Δx[i])
        integral += y[i + 1] * prefac * (Δx[i] + Δx[i + 1])^2 / (Δx[i] * Δx[i + 1])
        integral += y[i + 2] * prefac * (2 - Δx[i] / Δx[i + 1])
    end

    if isodd(N)
        integral += y[end] * (2 * Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) /
             (6 * (Δx[end - 1] + Δx[end]))
        integral += y[end - 1] * (Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) / (6 * Δx[end - 1])
        integral -= y[end - 2] * Δx[end]^3 / (6 * Δx[end - 1] * (Δx[end - 1] + Δx[end]))
    end

    return integral
end

function integrate(x::AbstractVector, y::AbstractVector, dx::AbstractVector, ::Union{Real,AbstractVector}, ::QESimpson)
    length(y) >= 4 || error("y must contain at least four elements")

    n = length(x)  # Number of points (number of intervals + 1)
    ifirst = firstindex(x)
    ilast = lastindex(x)

    integral = sum((ifirst + 1):(ilast - 1)) do i
        2 * abs(mod(i, 2) - 2) * y[i] * dx[i]
    end
    if mod(n, 2) == 1
        integral = (integral + y[ifirst] * dx[ifirst] + y[ilast  ] * dx[ilast  ]) / 3
    else
        integral = (integral + y[ifirst] * dx[ifirst] + y[ilast-1] * dx[ilast-1]) / 3
    end
    return integral
end

function integrate(::AbstractVector, y::AbstractVector, ::AbstractVector, Δx::Real, ::AbinitCorrectedTrapezoid)
    n = length(y)  # Number of points (number of intervals + 1)
    ifirst = firstindex(y)
    ilast = lastindex(y)

    if n >= 10
        endpoint = (
            23.75(y[ifirst    ] + y[ilast    ]) +
            95.10(y[ifirst + 1] + y[ilast - 1]) +
            55.20(y[ifirst + 2] + y[ilast - 2]) +
            79.30(y[ifirst + 3] + y[ilast - 3]) +
            70.65(y[ifirst + 4] + y[ilast - 4])
        ) / 72
        return (sum(y[ifirst+5:ilast-5]) + endpoint) * Δx
    elseif n >= 8  # https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_3/8_rule
        endpoint = (
            17(y[ifirst    ] + y[ilast    ]) +
            59(y[ifirst + 1] + y[ilast - 1]) +
            43(y[ifirst + 2] + y[ilast - 2]) +
            49(y[ifirst + 3] + y[ilast - 3])
        ) / 48
        return n == 8 ? endpoint * Δx : (endpoint + y[ifirst + 4]) * Δx
    elseif n == 7  # slightly modified version of n >=8 case
        return (
            17(y[ifirst    ] + y[ilast    ]) +
            59(y[ifirst + 1] + y[ilast - 1]) +
            43(y[ifirst + 2] + y[ilast - 2]) +
            50(y[ifirst + 3]                )
        ) / 48 * Δx
    elseif n == 6
        return (
            17(y[ifirst    ] + y[ilast    ]) +
            59(y[ifirst + 1] + y[ilast - 1]) +
            44(y[ifirst + 2] + y[ilast - 2])
        ) / 48 * Δx
    elseif n == 5  # Simpson's 1/3 rule
        return (
             (y[ifirst    ] + y[ifirst + 4]) +
            4(y[ifirst + 1] + y[ifirst + 3]) +
            2(y[ifirst + 2]                )
        ) / 3 * Δx
    elseif n == 4
        return (
            3(y[ifirst   ] + y[ifirst + 1]) +
            9(y[ilast - 1] + y[ilast     ])
        ) / 8 * Δx
    elseif n == 3  # Simpson's 1/3 rule
        return (
             (y[ifirst    ] + y[ifirst + 2]) +
            4(y[ifirst + 1] + y[ifirst + 1])
        ) / 3 * Δx
    elseif n == 2  # Trapezoidal rule
        return (y[ifirst] + y[ilast]) / 2 * Δx
    elseif n == 1
        return y[ifirst] * Δx
    end
end
