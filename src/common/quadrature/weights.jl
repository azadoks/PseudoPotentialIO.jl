function integration_weights(x::AbstractVector, dx::AbstractVector,
                             Δx::Union{Real,AbstractVector}, method::QuadratureMethod)
    weights = similar(x)
    return integration_weights!(weights, x, dx, Δx, method)
end

function integration_weights!(weights::AbstractVector, ::AbstractVector,
                              ::AbstractVector, Δx::Real, ::Trapezoid)
    weights[begin] = Δx / 2
    weights[(begin + 1):(end - 1)] .= Δx
    weights[end] = Δx / 2

    return weights
end

function integration_weights!(weights::AbstractVector, x::AbstractVector,
                              ::AbstractVector, Δx::AbstractVector, ::Trapezoid)
    weights[begin] = Δx[begin] / 2
    for i in (firstindex(weights) + 1):(lastindex(weights) - 1)
        # Δx[i] + Δx[i-1] = (x[i+1] - x[i]) + (x[i] - x[i-1]) = x[i+i] - x[i-1]
        weights[i] = (x[i + 1] - x[i - 1]) / 2
    end
    weights[end] = Δx[end] / 2

    return weights
end

function integration_weights!(weights::AbstractVector, x::AbstractVector,
                              ::AbstractVector, Δx::Real, ::Simpson)
    N = length(x) - 1  # Number of intervals

    if !isodd(N)  # Standard Simpson's composite 1/3 rule
        weights[begin] = 1 / 3 * Δx
        weights[(begin + 1):2:(end - 1)] .= 4 * 1 / 3 * Δx
        weights[(begin + 2):2:(end - 1)] .= 2 * 1 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    else # * For functions which decay at the right bound, this is a better approximation
        # If the number of intervals is odd, apply Simpsons method to the first N-1 intervals
        # and the Trapezoidal method to the last interval.
        weights[begin] = 1 / 3 * Δx
        weights[(begin+1):2:(end-2)] .= 4 * 1 / 3 * Δx
        weights[(begin+2):2:(end-2)] .= 2 * 1 / 3 * Δx
        weights[end-1] = 5 / 6 * Δx
        weights[end] = 1 / 2 * Δx
    end
    # else
    #     # If the number of intervals is odd, apply Simpsons method to the last N-1 intervals
    #     # and the Trapezoidal method to the first interval.
    #     weights[begin] = 1 / 2 * Δx
    #     weights[begin+1] = 5 / 6 * Δx
    #     weights[(begin+2):2:(end-1)] .= 4 * 1 / 3 * Δx
    #     weights[(begin+3):2:(end-1)] .= 2 * 1 / 3 * Δx
    #     weights[end] = 1 / 3 * Δx
    # end
    # else
    #     # If the number of intervals is odd, average the results of applying Simpson's
    #     # composite 1/3 rule to the first N-1 and last N-1 intervals, applying the
    #     # Trapezoidal rule to the last / first interval.
    #     weights[begin] = 5 / 12 * Δx  # (1/3 + 1/2) / 2 = 5/12
    #     weights[begin + 1] = 13 / 12 * Δx  # (4/3 + 1/3 + 1/2) / 2 = 13/12
    #     weights[(begin + 2):(end - 2)] .= 1 * Δx  # (4/3 + 2/3) / 2 = 1
    #     weights[end - 1] = 13 / 12 * Δx  # ''
    #     weights[end] = 5 / 12 * Δx  # ''
    # end

    return weights
end

function integration_weights!(weights::AbstractVector, x::AbstractVector,
                              ::AbstractVector, Δx::AbstractVector, ::Simpson)
    N = length(x) - 1  # Number of intervals
    fill!(weights, 0)

    # Skip the last interval if the number of intervals is odd
    istop = isodd(N) ? lastindex(x) - 3 : lastindex(x) - 2

    for i in firstindex(x):2:istop
        prefac = (Δx[i] + Δx[i + 1]) / 6
        weights[i] += prefac * (2 - Δx[i + 1] / Δx[i])
        weights[i + 1] += prefac * (Δx[i] + Δx[i + 1])^2 / (Δx[i] * Δx[i + 1])
        weights[i + 2] += prefac * (2 - Δx[i] / Δx[i + 1])
    end

    if isodd(N)  # This handles the last interval when the number of intervals is odd
        weights[end] += (2 * Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) /
                        (6 * (Δx[end - 1] + Δx[end]))
        weights[end - 1] += (Δx[end]^2 + 3 * Δx[end] * Δx[end - 1]) / (6 * Δx[end - 1])
        weights[end - 2] -= Δx[end]^3 / (6 * Δx[end - 1] * (Δx[end - 1] + Δx[end]))
    end

    return weights
end

function integration_weights!(weights::AbstractVector, x::AbstractVector,
                              dx::AbstractVector, ::Union{Real,AbstractVector},
                              ::QESimpson)
    i = (firstindex(x) + 1):(lastindex(x) - 1)
    weights[i] .= 2 / 3 * abs.(mod.(i, 2) .- 2) .* dx[i]
    if mod(length(x), 2) == 1
        weights[begin] = 1 / 3 * dx[begin]
        weights[end] = 1 / 3 * dx[end]
    else
        weights[begin] = 1 / 3 * dx[begin]
        weights[end - 1] = 1 / 3 * dx[end - 1]
    end
    return weights
end

function integration_weights!(weights::AbstractVector, x::AbstractVector,
                              ::AbstractVector, Δx::Real, ::AbinitCorrectedTrapezoid)
    if length(x) >= 10
        weights[begin] = 23.75 / 72 * Δx
        weights[begin + 1] = 95.10 / 72 * Δx
        weights[begin + 2] = 55.20 / 72 * Δx
        weights[begin + 3] = 79.30 / 72 * Δx
        weights[begin + 4] = 70.65 / 72 * Δx
        weights[(begin + 5):(end - 5)] .= Δx
        weights[end - 4] = 70.65 / 72 * Δx
        weights[end - 3] = 79.30 / 72 * Δx
        weights[end - 2] = 55.20 / 72 * Δx
        weights[end - 1] = 95.10 / 72 * Δx
        weights[end] = 23.75 / 72 * Δx
    elseif length(x) == 9
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[begin + 4] = Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(x) == 8
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(x) == 7
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 50 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(x) == 6
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 44 / 48 * Δx
        weights[end - 2] = 44 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(x) == 5
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 4 / 3 * Δx
        weights[begin + 2] = 2 / 3 * Δx
        weights[end - 1] = 4 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(x) == 4
        weights[begin] = 3 / 8 * Δx
        weights[begin + 1] = 9 / 8 * Δx
        weights[end - 1] = 9 / 8 * Δx
        weights[end] = 3 / 8 * Δx
    elseif length(x) == 3
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 8 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(x) == 2
        weights[begin] = 1 / 2 * Δx
        weights[end] = 1 / 2 * Δx
    elseif length(x) == 1
        weights[begin] = Δx
    end
    return weights
end
