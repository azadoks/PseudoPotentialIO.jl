function integrate(x::AbstractVector, f, ::AbstractVector, Δx::Real, ::Trapezoid)
    length(x) >= 2 || error("x and y must contain at least two elements")

    return (sum(f, x) - (f(x[begin]) + f(x[end])) / 2) * Δx
end

function integrate(x::AbstractVector, f, ::AbstractVector, Δx::AbstractVector, ::Trapezoid)
    length(x) - 1 == length(Δx) || error("Δx must contain one fewer element than x and y")
    length(x) >= 2 || error("x and y must contain at least two elements")

    ifirst = firstindex(x)
    ilast = lastindex(x)

    integral = (f(x[begin]) + f(x[begin + 1])) * Δx[begin]
    for i in (ifirst + 1):(ilast - 1)
        integral += (f(x[i]) + f(x[i + 1])) * Δx[i]
    end
    return integral / 2
end

function integrate(x::AbstractVector, f, ::AbstractVector, Δx::Real, ::Simpson)
    length(x) >= 4 || error("x and y must contain at least four elements")

    N = length(x) - 1  # Number of intervals
    if !isodd(N)
        integral = f(x[begin])
        integral += 4 * sum(f, x[(begin + 1):2:(end - 1)])
        integral += 2 * sum(f, x[(begin + 2):2:(end - 1)])
        integral += f(x[end])
        return integral / 3 * Δx
    else
        integral = 5 / 4 * f(x[begin])
        integral += 13 / 4 * f(x[begin + 1])
        integral += 3 * sum(f, x[(begin + 2):(end - 2)])
        integral += 13 / 4 * f(x[end - 1])
        integral += 5 / 4 * f(x[end - 1])
        return integral / 3 * Δx
    end
end

function integrate(x::AbstractVector, f, ::AbstractVector, Δx::AbstractVector, ::Simpson)
    length(x) - 1 == length(Δx) || error("Δx must contain one fewer element than x and y")
    length(x) >= 4 || error("x and y must contain at least four elements")

    N = length(x) - 1  # Number of intervals
    istop = isodd(N) ? lastindex(x) - 3 : lastindex(x) - 2

    integral = 0
    for i in firstindex(x):2:istop
        prefac = (Δx[i] + Δx[i + 1]) / 6
        integral += f(x[i]) * prefac * (2 - Δx[i + 1] / Δx[i])
        integral += f(x[i + 1]) * prefac * (Δx[i] + Δx[i + 1])^2 / (Δx[i] * Δx[i + 1])
        integral += f(x[i + 2]) * prefac * (2 - Δx[i] / Δx[i + 1])
    end

    if isodd(N)
        integral += f(x[end]) * (2 * Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) /
             (6 * (Δx[end - 1] + Δx[end]))
        integral += f(x[end - 1]) * (Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) / (6 * Δx[end - 1])
        integral -= f(x[end - 2]) * Δx[end]^3 / (6 * Δx[end - 1] * (Δx[end - 1] + Δx[end]))
    end

    return integral
end

function integrate(x::AbstractVector, f, dx::AbstractVector, ::Union{Real,AbstractVector}, ::QESimpson)
    length(x) >= 4 || error("x, y, and dx must contain at least four elements")

    n = length(x)  # Number of points (number of intervals + 1)

    integral = sum(firstindex(x)+1:lastindex(x)-1) do i
        2 * abs(mod(i, 2) - 2) * f(x[i]) * dx[i]
    end
    if mod(n, 2) == 1
        integral = (integral + f(x[begin]) * dx[begin] + f(x[end  ]) * dx[end  ]) / 3
    else
        integral = (integral + f(x[begin]) * dx[begin] + f(x[end-1]) * dx[end-1]) / 3
    end
    return integral
end

function integrate(x::AbstractVector, f, ::AbstractVector, Δx::Real, ::AbinitCorrectedTrapezoid)
    n = length(x)  # Number of points (number of intervals + 1)

    if n >= 10
        endpoint = (
            23.75(f(x[begin    ]) + f(x[end    ])) +
            95.10(f(x[begin + 1]) + f(x[end - 1])) +
            55.20(f(x[begin + 2]) + f(x[end - 2])) +
            79.30(f(x[begin + 3]) + f(x[end - 3])) +
            70.65(f(x[begin + 4]) + f(x[end - 4]))
        ) / 72
        return (sum(f, x[begin+5:end-5]; init=0) + endpoint) * Δx
    elseif n >= 8  # https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_3/8_rule
        endpoint = (
            17(f(x[begin    ]) + f(x[end    ])) +
            59(f(x[begin + 1]) + f(x[end - 1])) +
            43(f(x[begin + 2]) + f(x[end - 2])) +
            49(f(x[begin + 3]) + f(x[end - 3]))
        ) / 48
        return n == 8 ? endpoint * Δx : (endpoint + f(x[begin + 4])) * Δx
    elseif n == 7  # slightly modified version of n >=8 case
        return (
            17(f(x[begin    ]) + f(x[end    ])) +
            59(f(x[begin + 1]) + f(x[end - 1])) +
            43(f(x[begin + 2]) + f(x[end - 2])) +
            50(f(x[begin + 3])                )
        ) / 48 * Δx
    elseif n == 6
        return (
            17(f(x[begin    ]) + f(x[end    ])) +
            59(f(x[begin + 1]) + f(x[end - 1])) +
            44(f(x[begin + 2]) + f(x[end - 2]))
        ) / 48 * Δx
    elseif n == 5  # Simpson's 1/3 rule
        return (
             (f(x[begin    ]) + f(x[begin + 4])) +
            4(f(x[begin + 1]) + f(x[begin + 3])) +
            2(f(x[begin + 2])                  )
        ) / 3 * Δx
    elseif n == 4
        return (
            3(f(x[begin   ]) + f(x[begin + 1])) +
            9(f(x[end - 1 ]) + f(x[end      ]))
        ) / 8 * Δx
    elseif n == 3  # Simpson's 1/3 rule
        return (
             (f(x[begin    ]) + f(x[begin + 2])) +
            4(f(x[begin + 1]) + f(x[begin + 1]))
        ) / 3 * Δx
    elseif n == 2  # Trapezoidal rule
        return (f(x[begin]) + f(x[end])) / 2 * Δx
    elseif n == 1
        return f(x[begin]) * Δx
    end
end
