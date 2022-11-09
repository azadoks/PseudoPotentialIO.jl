struct NormConservingPsP{T} <: NumericPsP{T}
    "Total charge."
    Ztot::T
    "Valence charge."
    Zval::T
    "Maximum angular momentum."
    lmax::Int
    "Radial mesh."
    r::Vector{T}
    "Radial mesh spacing."
    dr::Union{T,Vector{T}}
    "Local part of the potential on the radial mesh."
    Vloc::Vector{T}
    "Nonlocal projectors β[l][n] on the radial mesh."
    β::OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}
    "Cutoff indices for nonlocal projectors"
    β_ircut::OffsetVector{Vector{Int},Vector{Vector{Int}}}
    "Projector coupling coefficients D[l][n,m]."
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions ϕ̃[l][n] on the radial mesh."
    ϕ̃::Union{Nothing,OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}}
    "Cutoff indices for nonlocal projectors"
    ϕ̃_ircut::Union{Nothing,OffsetVector{Vector{Int},Vector{Vector{Int}}}}
    "Model core charge density for non-linear core correction on the radial mesh."
    ρcore::Union{Nothing,Vector{T}}
    "Valence charge density for charge density initialization on the radial mesh."
    ρval::Union{Nothing,Vector{T}}
end

function NormConservingPsP(upf::UpfFile)
    if (upf.header.pseudo_type != "NC")
        error("Provided `UpfFile` is not a norm-conserving pseudo")
    end
    if (upf.header.relativistic == "full") | upf.header.has_so | (!isnothing(upf.spin_orb))
        error("Fully relativistic pseudos are not supported")
    end
    return _construct_nc_internal(upf)
end

function _construct_nc_internal(upf::UpfFile)
    Ztot = PeriodicTable.elements[Symbol(upf.header.element)].number
    Zval = upf.header.z_valence
    lmax = upf.header.l_max
    r = upf.mesh.r
    Vloc = upf.local_ ./ 2  # Ry -> Ha
    ρcore = isnothing(upf.nlcc) ? nothing : upf.nlcc
    # Indices in upf.nonlocal.betas for projectors at each angular momentum
    iβ_upf = OffsetArray(map(l -> filter(i -> upf.nonlocal.betas[i].angular_momentum == l,
                                         eachindex(upf.nonlocal.betas)), 0:lmax), 0:lmax)
    # Number of projectors at each angular momentum
    nβ = OffsetArray(length.(iβ_upf), 0:lmax)
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the 
    # sub-arrays D[l][n,m] can be extracted
    cum_nβ = [0, cumsum(nβ)...]
    # Extract the blocks from `upf.nonlocal.dij` corresponding to each angular momentum
    D = OffsetVector(map(i -> collect(upf.nonlocal.dij[(cum_nβ[i] + 1):cum_nβ[i + 1],
                                                       (cum_nβ[i] + 1):cum_nβ[i + 1]]),
                         1:(length(cum_nβ) - 1)), 0:lmax) ./ 2  # Ry -> Ha

    mesh_type, _, _ = guess_mesh_type(r, upf.mesh.rab)
    mesh_type == "unknown" && error("Unknown mesh type")

    β_ircut = OffsetVector(map(l -> map(i -> length(upf.nonlocal.betas[i].beta), iβ_upf[l]),
                               0:lmax), 0:lmax)
    if mesh_type in ("log_1", "log_2")
        dr = upf.mesh.rab
        β = OffsetVector(map(l -> map(i -> _standardize_β_ϕ̃(upf.nonlocal.betas[i].beta, r),
                                      iβ_upf[l]), 0:lmax), 0:lmax)
        ρval = _standardize_ρval(upf.rhoatom, r)
    else  # mesh_type == "linear"
        dr = upf.mesh.rab[1]
        β = OffsetVector(map(l -> map(i -> _extrapolate_standardize_β_ϕ̃(upf.nonlocal.betas[i].beta,
                                                                         r), iβ_upf[l]),
                             0:lmax), 0:lmax)
        ρval = _extrapolate_standardize_ρval(upf.rhoatom, r)
    end

    if !isnothing(upf.pswfc)
        # Collect the indices in upf.nonlocal.betas for projectors at each angular momentum
        iχ_upf = OffsetArray(map(l -> filter(i -> upf.pswfc[i].l == l, eachindex(upf.pswfc)),
                                 0:lmax), 0:lmax)
        ϕ̃_ircut = OffsetVector(map(l -> map(i -> length(upf.pswfc[i].chi), iχ_upf[l]),
                                    0:lmax), 0:lmax)
        if mesh_type in ("log_1", "log_2")
            ϕ̃ = OffsetVector(map(l -> map(i -> _standardize_β_ϕ̃(upf.pswfc[i].chi, r),
                                           iχ_upf[l]), 0:lmax), 0:lmax)
        else  # mesh_type == "linear"
            ϕ̃ = OffsetVector(map(l -> map(i -> _extrapolate_standardize_β_ϕ̃(upf.pswfc[i].chi,
                                                                              r), iχ_upf[l]),
                                  0:lmax), 0:lmax)
        end
    else
        ϕ̃_ircut = nothing
        ϕ̃ = nothing
    end

    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, β_ircut, D, ϕ̃,
                                      ϕ̃_ircut, ρcore, ρval)
end

function _standardize_β_ϕ̃(f_upf::Vector{Float64}, r::Vector{Float64})::Vector{Float64}
    return @views f_upf ./ r[1:length(f_upf)]
end

function _standardize_ρval(ρval_upf::Vector{Float64}, r::Vector{Float64})::Vector{Float64}
    return @. ρval_upf / (4Float64(π) * r^2)
end

function _extrapolate_standardize_β_ϕ̃(f_upf::Vector{Float64},
                                       r::Vector{Float64})::Vector{Float64}
    ircut = length(f_upf)
    f = zeros(Float64, ircut)
    f[2:ircut] = @views f_upf[2:ircut] ./ r[2:ircut]
    f[1] = fit(r[2:5], f[2:5])(r[1])
    return f
end

function _extrapolate_standardize_ρval(ρval_upf::Vector{Float64},
                                       r::Vector{Float64})::Vector{Float64}
    ρval = zeros(Float64, length(ρval_upf))
    ρval[2:end] = @views ρval_upf[2:end] ./ (4Float64(π) * r[2:end] .^ 2)
    ρval[1] = fit(r[2:5], ρval[2:5])(r[1])
    return ρval
end

function NormConservingPsP(psp8::Psp8File)
    if psp8.header.extension_switch in (2, 3)
        error("Fully relativistic pseudos are not supported")
    end

    Ztot = psp8.header.zatom
    Zval = psp8.header.zion
    lmax = psp8.header.lmax
    r = psp8.rgrid
    dr = mean(diff(r))
    Vloc = psp8.v_local
    #TODO find β ircut, cut off β
    β = OffsetVector(map(l -> map(i -> psp8.projectors[l + 1][i],
                                  eachindex(psp8.projectors[l + 1])), 0:lmax), 0:lmax)
    β_ircut = OffsetVector(map(l -> map(i -> length(β[l][i]), eachindex(β[l])), 0:lmax),
                           0:lmax)
    D = OffsetVector(map(l -> diagm(psp8.ekb[l + 1]), 0:lmax), 0:lmax)
    ϕ̃ = nothing
    ϕ̃_ircut = nothing
    ρcore = psp8.rhoc
    ρval = nothing
    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, β_ircut, D, ϕ̃,
                                      ϕ̃_ircut, ρcore, ρval)
end

formalism(::NormConservingPsP)::Symbol = :norm_conserving
