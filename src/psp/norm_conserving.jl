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
    # There are two possible units schemes for the projectors and coupling coefficients:
    # β [Ry Bohr^{-1/2}]  D [Ry^{-1}]
    # β [Bohr^{-1/2}]     D [Ry]
    # The quantity that's used in calculations is β D β, so the units don't practically
    # matter. However, HGH pseudos in UPF format use the first units, so we assume them
    # to facilitate comparison of the intermediate quantities with analytical HGH.

    Ztot = PeriodicTable.elements[Symbol(upf.header.element)].number
    Zval = upf.header.z_valence
    lmax = upf.header.l_max
    r = upf.mesh.r
    Vloc = upf.local_ ./ 2  # Ry -> Ha
    ρcore = isnothing(upf.nlcc) ? nothing : upf.nlcc
    
    # Indices in upf.nonlocal.betas for projectors at each angular momentum
    iβ_upf = map(0:lmax) do l
        β_upf = upf.nonlocal.betas
        filter(i -> β_upf[i].angular_momentum == l, eachindex(β_upf))
    end
    iβ_upf = OffsetVector(iβ_upf, 0:lmax)
    
    # Number of projectors at each angular momentum
    nβ = OffsetArray(length.(iβ_upf), 0:lmax)
    
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the 
    # sub-arrays D[l][n,m] can be extracted
    cumul_nβ = [0, cumsum(nβ)...]
    
    # Extract the blocks from `upf.nonlocal.dij` corresponding to each angular momentum
    D = map(1:length(cumul_nβ) - 1) do i
        collect(upf.nonlocal.dij[(cumul_nβ[i] + 1):cumul_nβ[i + 1],
                                 (cumul_nβ[i] + 1):cumul_nβ[i + 1]])
    end
    D = OffsetVector(D, 0:lmax) .* 2  # 1/Ry -> 1/Ha

    # Guess the mesh type to know how to extrapolate to r=0 when converting the projectors
    # from the file (rβ) to the standard quantity (β)
    mesh_type, _, _ = guess_mesh_type(r, upf.mesh.rab)
    mesh_type == "unknown" && error("Unknown mesh type")
    dr = mesh_type == "linear" ? upf.mesh.rab[1] : upf.mesh.rab

    # Find the cutoff radius index for each projector
    β_ircut = map(0:lmax) do l
        map(iβ_upf[l]) do i
            length(upf.nonlocal.betas[i].beta)
        end
    end
    β_ircut = OffsetVector(β_ircut, 0:lmax)

    # UPFs store the projectors multiplied by the radial grid. For compatability with the
    # PseudoPotentialIO interface, we need the pure projectors. To get them, we divide
    # rβ by r where r > 0 and use that data to extrapolate to r = 0 using a quadratic
    # polynomial.
    β = map(0:lmax) do l
        map(iβ_upf[l]) do i
            if upf.mesh.r[1] == 0
                _extrapolate_standardize_β_ϕ̃(upf.nonlocal.betas[i].beta, r, mesh_type)
            else
                _standardize_β_ϕ̃(upf.nonlocal.betas[i].beta, r)
            end
        end
    end
    β = OffsetVector(β, 0:lmax) ./ 2  # Ry -> Ha

    # UPFs store the pseudo-atomic valence charge density with a prefactor of 4π multiplied
    # by the square of the radial grid. We do the extrapolation like for the projectors.
    if upf.mesh.r[1] == 0
        ρval = _extrapolate_standardize_ρval(upf.rhoatom, r, mesh_type)
    else
        ρval = _standardize_β_ϕ̃(upf.rhoatom, r)
    end

    # The pseudo-atomic wavefunctions need the same extrapolation treatment as the
    # projectors.
    if !isnothing(upf.pswfc)
        # Collect the indices in upf.nonlocal.betas for projectors at each angular momentum
        iχ_upf = map(0:lmax) do l
            filter(i -> upf.pswfc[i].l == l, eachindex(upf.pswfc))
        end
        iχ_upf = OffsetVector(iχ_upf, 0:lmax)
        
        ϕ̃_ircut = map(0:lmax) do l
            map(iχ_upf[l]) do i
                length(upf.pswfc[i].chi)
            end
        end
        ϕ̃_ircut = OffsetVector(ϕ̃_ircut, 0:lmax)

        ϕ̃ = map(0:lmax) do l
            map(iχ_upf[l]) do i
                if upf.mesh.r[1] == 0
                    _extrapolate_standardize_β_ϕ̃(upf.pswfc[i].chi, r, mesh_type)
                else
                    _standardize_β_ϕ̃(upf.pswfc[i].chi, r)
                end
            end
        end
        ϕ̃ = OffsetVector(ϕ̃, 0:lmax)
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
                                      r::Vector{Float64}, mesh_type::String)::Vector{Float64}
    ircut = length(f_upf)
    f = zeros(Float64, ircut)
    if mesh_type in ("log_1", "log_2")
        f[3:ircut] = @views f_upf[3:ircut] ./ r[3:ircut]
        poly = fit(r[3:6], f[3:6], 2)
        f[1:2] = poly.(r[1:2])
    else
        f[2:ircut] = @views f_upf[2:ircut] ./ r[2:ircut]
        f[1] = fit(r[2:5], f[2:5], 2)(r[1])
    end
    f
    return f
end

function _extrapolate_standardize_ρval(ρval_upf::Vector{Float64},
                                       r::Vector{Float64}, mesh_type::String)::Vector{Float64}
    ρval = zeros(Float64, length(ρval_upf))
    if mesh_type in ("log_1", "log_2")
        ρval[3:end] = @views ρval_upf[3:end] ./ (4Float64(π) * r[3:end] .^ 2)
        poly = fit(r[3:6], ρval[3:6], 2)
        ρval[1:2] = poly.(r[1:2])
    else
        ρval[2:end] = @views ρval_upf[2:end] ./ (4Float64(π) * r[2:end] .^ 2)
        ρval[1] = fit(r[2:5], ρval[2:5], 2)(r[1])
    end
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
