struct UltrasoftPsP{T} <: NumericPsP{T}
    "Total charge"
    Ztot::T
    "Valence charge"
    Zval::T
    "Maximum angular momentum"
    lmax::Int
    "Radial mesh"
    r::Vector{T}
    "Radial mesh spacing"
    dr::Union{T,Vector{T}}
    "Local part of the potential on the radial mesh"
    Vloc::Vector{T}
    "Nonlocal projectors β[l][n] on the radial mesh"
    β::OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
    "Cutoff indices for nonlocal projectors"
    β_ircut::OffsetVector{Vector{Int},Vector{Vector{Int}}}
    "Projector coupling coefficients D[l][n,m]"
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions ϕ̃[l][n] on the radial mesh."
    ϕ̃::OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
    "Cutoff indices for nonlocal projectors"
    ϕ̃_ircut::OffsetVector{Vector{Int},Vector{Vector{Int}}}
    "Augmentation charge density functions Q[l][n,m] on the radial mesh"
    Q::OffsetVector{Matrix{Vector{T}},Vector{Matrix{Vector{T}}}}
    "Augmentation charges q[l][n,m]"
    q::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Model core charge density for non-linear core correction on the radial mesh"
    ρcore::Union{Nothing,Vector{T}}
    "Valence charge density for charge density initialization on the radial mesh"
    ρval::Union{Nothing,Vector{T}}
end

function UltrasoftPsP(upf::UpfFile)
    if !in(upf.header.pseudo_type, ("US", "USPP"))
        error("Provided `UpfFile` is not an ultrasoft pseudo")
    end
    if (upf.header.relativistic == "full") | upf.header.has_so | (!isnothing(upf.spin_orb))
        error("Fully relativistic pseudos are not supported")
    end
    return _upf_construct_us_internal(upf)
end

function _upf_construct_us_internal(upf::UpfFile)
    nc = _upf_construct_nc_internal(upf)
    # Number of projectors at each angular momentum
    nβ = OffsetVector(length.(nc.β), 0:(nc.lmax))
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the 
    # sub-arrays q[l][n,m] can be extracted
    cum_nβ = [0, cumsum(nβ)...]
    q_upf = upf.nonlocal.augmentation.q
    q = OffsetVector(map(i -> collect(q_upf[(cum_nβ[i] + 1):cum_nβ[i + 1],
                                            (cum_nβ[i] + 1):cum_nβ[i + 1]]),
                         1:(length(cum_nβ) - 1)), 0:(nc.lmax))
    # Wrangle the agumentation functions. For UPF v2.0.1, they are given with indeices
    # i ∈ 1:nβ, j ∈ 1:nβ, l ∈ 0:2lmax. In the old format, they are given only at i and j
    # with additional coefficients for a polynomial expansion within a cutoff radius
    # These are used to reconstruct the full Qijl. 
    if upf.nonlocal.augmentation.q_with_l
        Q = OffsetVector([Matrix{Vector{Float64}}(undef, upf.header.number_of_proj,
                                                  upf.header.number_of_proj)
                          for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
        for l in 0:(2upf.header.l_max), i in 1:(upf.header.number_of_proj),
            j in 1:(upf.header.number_of_proj)

            Q[l][i, j] = zeros(length(upf.mesh.r))
        end
        for l in 0:(2upf.header.l_max)
            Q_upf_l = filter(qijl -> qijl.angular_momentum == l,
                             upf.nonlocal.augmentation.qijls)
            for Q_upf in Q_upf_l
                Q[l][Q_upf.first_index, Q_upf.second_index] = Q_upf.qijl
                Q[l][Q_upf.second_index, Q_upf.first_index] = Q_upf.qijl
            end
        end
    elseif upf.nonlocal.augmentation.nqf > 0
        #TODO check correctness
        r = upf.mesh.r
        nqf = upf.nonlocal.augmentation.nqf
        nqlc = 2upf.header.l_max + 1

        Q = OffsetVector([Matrix{Vector{Float64}}(undef, upf.header.number_of_proj,
                                                  upf.header.number_of_proj)
                          for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
        for l in 0:(2upf.header.l_max), i in 1:(upf.header.number_of_proj),
            j in 1:(upf.header.number_of_proj)

            Q[l][i, j] = zeros(length(upf.mesh.r))
        end
        for (Q_upf, Qfcoef_upf) in
            zip(upf.nonlocal.augmentation.qijs, upf.nonlocal.augmentation.qfcoefs)
            qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
            rinner = upf.nonlocal.augmentation.rinner

            i = Q_upf.first_index
            j = Q_upf.second_index

            li = upf.nonlocal.betas[i].angular_momentum
            lj = upf.nonlocal.betas[j].angular_momentum

            for l in abs(li - lj):2:(li + lj)
                #* Reference implementation
                # qij = copy(Q_upf.qij)
                # for ir in eachindex(r) 
                #     if r[ir] < rinner[l + 1]
                #         qij[ir] = qfcoef[1, l + 1]
                #         for n in 2:nqf
                #             qij[ir] += qfcoef[n, l + 1] * r[ir]^(2n)
                #         end
                #         qij[ir] *= r[ir]^(l + 2)
                #     end
                # end

                #TODO not sure why this isn't equivalent
                # poly = Polynomial(qfcoef[:, l + 1])
                # qij[1:ircut] = r[1:ircut].^(l + 2) .* poly.(r[1:ircut].^2)

                qij = copy(Q_upf.qij)
                ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1

                qij[1:ircut] .= qfcoef[1, l + 1]
                for n in 2:nqf
                    qij[1:ircut] .+= qfcoef[n, l + 1] .* r[1:ircut] .^ (2n)
                end
                qij[1:ircut] .*= r[1:ircut] .^ (l + 2)

                Q[l][Q_upf.first_index, Q_upf.second_index] = qij
                Q[l][Q_upf.second_index, Q_upf.first_index] = qij
            end
        end
    else
        error("q_with_l = false and nqf == 0, unsure what to do...")
    end
    return UltrasoftPsP{Float64}(nc.Ztot, nc.Zval, nc.lmax, nc.r, nc.dr, nc.Vloc, nc.β,
                                 nc.β_ircut, nc.D, nc.ϕ̃, nc.ϕ̃_ircut, Q, q, nc.ρcore,
                                 nc.ρval)
end

is_ultrasoft(::UltrasoftPsP)::Bool = true

#TODO test the augmentation functions

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.q[l]
end

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int, n::Int)::T where {T<:Real}
    return psp.q[l][n, n]
end

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int, n::Int,
                               m::Int)::T where {T<:Real}
    return psp.q[l][n, m]
end

function augmentation_real(psp::UltrasoftPsP, l::Int, n::Int, m::Int,
                           r::T)::T where {T<:Real}
    return interpolate((psp.r,), psp.Q[l][n, m], (Gridded(Linear()),))(r)
end

function augmentation_real(psp::UltrasoftPsP, l::Int, n::Int, m::Int,
                           R::AbstractVector{T})::T where {T<:Real}
    return augmentation_real(psp, l, n, m, norm(R))
end

function augmentation_fourier(psp::UltrasoftPsP, l::Int, n::Int, m::Int,
                              q::T)::T where {T<:Real}
    f = @. psp.r^2 * fast_sphericalbesselj0(q * r) * psp.Q[l][n, m]
    return 4π * trapezoid(f, psp.dr)
end

function augmentation_fourier(psp::UltrasoftPsP, l::Int, n::Int, m::Int,
                              K::AbstractVector{T})::T where {T<:Real}
    return augmentation_fourier(psp, l, n, m, norm(K))
end
