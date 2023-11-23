function integration_weights(mesh::RadialMesh, method::QuadratureMethod)
    weights = similar(mesh)
    return integration_weights!(weights, mesh, method)
end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::Trapezoid)
    Δx = mesh.a

    weights[begin] = Δx / 2
    weights[(begin + 1):(end - 1)] .= Δx
    weights[end] = Δx / 2

    return weights
end

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::Trapezoid)
    weights[begin] = diff(mesh, 1) / 2
    for i in (firstindex(weights) + 1):(lastindex(weights) - 1)
        # Δx[i] + Δx[i-1] = (x[i+1] - x[i]) + (x[i] - x[i-1])
        #                 = x[i+i] - x[i-1]
        weights[i] = (mesh[i + 1] - mesh[i - 1]) / 2
    end
    weights[end] = diff(mesh, mesh.n - 1) / 2

    return weights
end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::Simpson)
    Δx = mesh.a
    N = length(weights) - 1  # Number of intervals

    if !isodd(N)  # Standard Simpson's composite 1/3 rule
        weights[begin] = 1 / 3 * Δx
        weights[(begin + 1):2:(end - 1)] .= 4 * 1 / 3 * Δx
        weights[(begin + 2):2:(end - 1)] .= 2 * 1 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    else # * For functions which decay to zero as r -> ∞, this is a better approximation
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

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::Simpson)
    N = length(weights) - 1  # Number of intervals
    fill!(weights, 0)

    # Skip the last interval if the number of intervals is odd
    istop = isodd(N) ? lastindex(weights) - 3 : lastindex(weights) - 2

    for i in firstindex(weights):2:istop
        Δx_0 = diff(mesh, i)
        Δx_1 = diff(mesh, i + 1)
        prefac = (Δx_0 + Δx_1) / 6
        weights[i] += prefac * (2 - Δx_1 / Δx_0)
        weights[i + 1] += prefac * (Δx_0 + Δx_1)^2 / (Δx_0 * Δx_1)
        weights[i + 2] += prefac * (2 - Δx_0 / Δx_1)
    end

    if isodd(N)  # This handles the last interval when the number of intervals is odd
        Δx_n = diff(mesh, mesh.n - 1)
        Δx_nm1 = diff(mesh, mesh.n - 2)
        weights[end] += (2 * Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * (Δx_nm1 + Δx_n))
        weights[end - 1] += (Δx_n^2 + 3 * Δx_n * Δx_nm1) / (6 * Δx_nm1)
        weights[end - 2] -= Δx_n^3 / (6 * Δx_nm1 * (Δx_nm1 + Δx_n))
    end

    return weights
end

function integration_weights!(weights::AbstractVector, mesh::RadialMesh, ::QESimpson)
    i = (firstindex(weights) + 1):(lastindex(weights) - 1)
    weights[i] .= 2 / 3 * abs.(mod.(i, 2) .- 2) .* deriv.(Ref(mesh), i)
    if mod(length(weights), 2) == 1
        weights[begin] = 1 / 3 * deriv(mesh, 1)
        weights[end] = 1 / 3 * deriv(mesh, mesh.n - 1)
    else
        weights[begin] = 1 / 3 * deriv(mesh, 1)
        weights[end - 1] = 1 / 3 * deriv(mesh, mesh.n - 2)
    end
    return weights
end

function integration_weights!(weights::AbstractVector, mesh::UniformMesh, ::AbinitCorrectedTrapezoid)
    Δx = mesh.a

    if length(weights) >= 10
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
    elseif length(weights) == 9
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[begin + 4] = Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 8
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 49 / 48 * Δx
        weights[end - 3] = 49 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 7
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 43 / 48 * Δx
        weights[begin + 3] = 50 / 48 * Δx
        weights[end - 2] = 43 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 6
        weights[begin] = 17 / 48 * Δx
        weights[begin + 1] = 59 / 48 * Δx
        weights[begin + 2] = 44 / 48 * Δx
        weights[end - 2] = 44 / 48 * Δx
        weights[end - 1] = 59 / 48 * Δx
        weights[end] = 17 / 48 * Δx
    elseif length(weights) == 5
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 4 / 3 * Δx
        weights[begin + 2] = 2 / 3 * Δx
        weights[end - 1] = 4 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(weights) == 4
        weights[begin] = 3 / 8 * Δx
        weights[begin + 1] = 9 / 8 * Δx
        weights[end - 1] = 9 / 8 * Δx
        weights[end] = 3 / 8 * Δx
    elseif length(weights) == 3
        weights[begin] = 1 / 3 * Δx
        weights[begin + 1] = 8 / 3 * Δx
        weights[end] = 1 / 3 * Δx
    elseif length(weights) == 2
        weights[begin] = 1 / 2 * Δx
        weights[end] = 1 / 2 * Δx
    elseif length(weights) == 1
        weights[begin] = Δx
    end
    return weights
end
