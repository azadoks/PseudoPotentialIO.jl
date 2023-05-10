"""
Type representing a numeric ultrasoft pseudopotential.
"""
struct UltrasoftPsP{T} <: NumericPsP{T}
    "Identifier"
    identifier::String
    "SHA1 Checksum"
    checksum::Union{Nothing,Vector{UInt8}}
    "Total charge"
    Zatom::Int
    "Valence charge"
    Zval::Int
    "Maximum angular momentum"
    lmax::Int
    "Radial mesh"
    r::Vector{T}
    "Radial mesh spacing"
    dr::Vector{T}
    Δr::Union{T,Vector{T}}
    "Local part of the potential on the radial mesh"
    Vloc::Vector{T}
    "Nonlocal projectors β[l][n] on the radial mesh"
    β::OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
    "Projector coupling coefficients D[l][n,m]"
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions χ[l][n] on the radial mesh."
    χ::OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
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

function _upf_construct_augmentation_q_with_l(upf::UpfFile)
    Q = OffsetVector([Matrix{Vector{Float64}}(undef, upf.header.number_of_proj,
                                              upf.header.number_of_proj)
                      for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
    for l in 0:(2upf.header.l_max)
        # Fill Q with zero vectors
        for i in 1:(upf.header.number_of_proj), j in 1:(upf.header.number_of_proj)
            Q[l][i, j] = zeros(length(upf.mesh.r))
        end
        # Replace the zero vectors with data from the UPF where present
        Q_upf_l = filter(qijl -> qijl.angular_momentum == l,
                         upf.nonlocal.augmentation.qijls)
        for Q_upf in Q_upf_l
            Q[l][Q_upf.first_index, Q_upf.second_index] = Q_upf.qijl
            Q[l][Q_upf.second_index, Q_upf.first_index] = Q_upf.qijl
        end
    end
    return Q
end

@views function _upf_construct_augmentation_qfcoef(upf::UpfFile)
    #TODO check correctness
    r = upf.mesh.r
    r2 = upf.mesh.r .^ 2
    nqf = upf.nonlocal.augmentation.nqf
    nqlc = 2upf.header.l_max + 1

    Q = OffsetVector([Matrix{Vector{Float64}}(undef, upf.header.number_of_proj,
                                              upf.header.number_of_proj)
                      for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
    for l in 0:(2upf.header.l_max), i in 1:(upf.header.number_of_proj),
        j in 1:(upf.header.number_of_proj)
        # Fill Q with zero vectors
        Q[l][i, j] = zeros(length(upf.mesh.r))
    end
    for (Q_upf, Qfcoef_upf) in
        zip(upf.nonlocal.augmentation.qijs, upf.nonlocal.augmentation.qfcoefs)
        # Replace the zero vectors with datat from the UPF where present
        # It's not worth the effort to make these into OffsetVectors zero-indexed for l.
        qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
        rinner = upf.nonlocal.augmentation.rinner

        i = Q_upf.first_index
        j = Q_upf.second_index

        li = upf.nonlocal.betas[i].angular_momentum
        lj = upf.nonlocal.betas[j].angular_momentum

        for l in abs(li - lj):2:(li + lj)
            qij = copy(Q_upf.qij)
            ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1
            poly = Polynomial(qfcoef[:, l + 1])
            qij[1:ircut] = r[1:ircut] .^ (l + 2) .* poly.(r2[1:ircut])

            Q[l][Q_upf.first_index, Q_upf.second_index] = qij
            Q[l][Q_upf.second_index, Q_upf.first_index] = qij
        end
    end
    return Q
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
        Q = _upf_construct_augmentation_q_with_l(upf)
    elseif upf.nonlocal.augmentation.nqf > 0
        Q = _upf_construct_augmentation_qfcoef(upf)
    else
        error("q_with_l == false and nqf == 0, unsure what to do...")
    end
    return UltrasoftPsP{Float64}(nc.identifier, nc.checksum, nc.Zatom, nc.Zval, nc.lmax,
                                 nc.r, nc.dr, nc.Δr, nc.Vloc, nc.β, nc.D, nc.χ, Q, q,
                                 nc.ρcore, nc.ρval)
end

is_norm_conserving(::UltrasoftPsP)::Bool = false
is_ultrasoft(::UltrasoftPsP)::Bool = true
is_paw(::UltrasoftPsP)::Bool = false

#TODO test the augmentation functions
has_quantity(::AugmentationCoupling, psp::UltrasoftPsP) = true
get_quantity(::AugmentationCoupling, psp::UltrasoftPsP, l) = psp.q[l]
get_quantity(::AugmentationCoupling, psp::UltrasoftPsP, l, n) = psp.q[l][n, n]
get_quantity(::AugmentationCoupling, psp::UltrasoftPsP, l, n, m) = psp.q[l][n, m]

function psp_quantity_evaluator(::RealSpace, q::AugmentationFunction, psp::UltrasoftPsP, l,
                                n, m)
    return build_interpolator_real(psp.Q[l][n, m], psp.r)
end

function psp_quantity_evaluator(::FourierSpace, q::AugmentationFunction, psp::UltrasoftPsP,
                                l, n, m)
    return hankel_transform(psp.Q[l][n, m], l, psp.r, psp.dr)
end
