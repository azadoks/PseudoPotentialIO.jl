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
    Q::OffsetVector{Matrix{Vector{T}}, Vector{Matrix{Vector{T}}}}
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
    return _construct_us_internal(upf)
end

function _construct_us_internal(upf::UpfFile)
    nc = _construct_nc_internal(upf)
    # Number of projectors at each angular momentum
    nβ = OffsetVector(length.(nc.β), 0:nc.lmax)
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the 
    # sub-arrays q[l][n,m] can be extracted
    cum_nβ = [0, cumsum(nβ)...]
    q_upf = upf.nonlocal.augmentation.q
    q = OffsetVector(map(i -> collect(q_upf[(cum_nβ[i] + 1):cum_nβ[i + 1],
                                            (cum_nβ[i] + 1):cum_nβ[i + 1]]),
                         1:(length(cum_nβ) - 1)), 0:nc.lmax)
    #TODO: reconstruct Q(r) for r < rinner for UPF v1.old
    if upf.nonlocal.augmentation.q_with_l
        Q = OffsetVector([Matrix{Vector{Float64}}(undef, sum(nβ), sum(nβ)) for l in 0:nc.lmax], 0:nc.lmax)
        for l in 0:nc.lmax
            Q_upf_l = filter(qijl -> qijl.angular_momentum == l, upf.nonlocal.augmentation.qijls)
            for Q_upf in Q_upf_l
                Q[l][Q_upf.first_index, Q_upf.second_index] = Q_upf.qijl
            end
        end
    else
        Q = OffsetVector([Matrix{Vector{Float64}}(undef, sum(nβ), sum(nβ)) for l in 0:0], 0:0)
        for Q_upf in upf.nonlocal.augmentation.qijs
            Q[0][Q_upf.first_index, Q_upf.second_index] = Q_upf.qij
        end
    end
    return UltrasoftPsP{Float64}(nc.Ztot, nc.Zval, nc.lmax, nc.r, nc.dr, nc.Vloc, nc.β,
                                 nc.β_ircut, nc.D, nc.ϕ̃, nc.ϕ̃_ircut, Q, q, nc.ρcore, nc.ρval)
end

formalism(::UltrasoftPsP)::Symbol = :ultrasoft

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int)::Matrix{T} where {T<:Real}
    return psp.q[l]
end

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int, n::Int)::T where {T<:Real}
    return psp.q[l][n, n]
end

function augmentation_coupling(psp::UltrasoftPsP{T}, l::Int, n::Int, m::Int)::T where {T<:Real}
    return psp.q[l][n, m]
end

function augmentation_fourier(psp::UltrasoftPsP, l::Int, n::Int, q::T)::T where {T<:Real}
    f = @. psp.r^2 * fast_sphericalbesselj0(q * r) * psp.Q[l][n]
    return 4π * trapezoid(f, psp.dr)
end
