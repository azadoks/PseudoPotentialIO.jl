@doc raw"""
Type representing a numeric norm-conserving pseudopotential.
"""
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
    "Projector coupling coefficients D[l][n,m]."
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions ϕ̃[l][n] on the radial mesh."
    ϕ̃::Union{Nothing,OffsetVector{Vector{Vector{T}},Vector{Vector{Vector{T}}}}}
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
    return _upf_construct_nc_internal(upf)
end

function _upf_construct_nc_internal(upf::UpfFile)
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

    # Guess the mesh type to choose scalar or vector `dr`
    mesh_type, _, _ = guess_mesh_type(r, upf.mesh.rab)
    mesh_type == "unknown" && error("Unknown mesh type")
    dr = mesh_type == "linear" ? upf.mesh.rab[1] : upf.mesh.rab

    # Indices in upf.nonlocal.betas for projectors at each angular momentum
    iβ_upf = map(0:lmax) do l
        β_upf = upf.nonlocal.betas
        return filter(i -> β_upf[i].angular_momentum == l, eachindex(β_upf))
    end
    iβ_upf = OffsetVector(iβ_upf, 0:lmax)

    # Number of projectors at each angular momentum
    nβ = OffsetArray(length.(iβ_upf), 0:lmax)

    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the 
    # sub-arrays D[l][n,m] can be extracted
    cumul_nβ = [0, cumsum(nβ)...]

    # Extract the blocks from `upf.nonlocal.dij` corresponding to each angular momentum
    D = map(1:(length(cumul_nβ) - 1)) do i
        return collect(upf.nonlocal.dij[(cumul_nβ[i] + 1):cumul_nβ[i + 1],
                                        (cumul_nβ[i] + 1):cumul_nβ[i + 1]])
    end
    D = OffsetVector(D, 0:lmax) .* 2  # 1/Ry -> 1/Ha

    # UPFs store the projectors multiplied by the radial grid, so we multiply again by the
    # grid for consistency.
    β = map(0:lmax) do l
        map(iβ_upf[l]) do n
            βln = upf.nonlocal.betas[n].beta
            βln = βln .* @view r[1:length(βln)]
        end
    end
    β = OffsetVector(β, 0:lmax) ./ 2  # Ry -> Ha

    # UPFs store the pseudo-atomic valence charge density with a prefactor of 4π r^2.
    # For consistency, we remove the 4π prefactor.
    ρval = upf.rhoatom ./ 4π

    if !isnothing(upf.pswfc)
        # Collect the indices in upf.pswfc for projectors at each angular momentum
        iχ_upf = map(0:lmax) do l
            return filter(i -> upf.pswfc[i].l == l, eachindex(upf.pswfc))
        end
        iχ_upf = OffsetVector(iχ_upf, 0:lmax)

        # UPFs store the wavefunctions multiplied by the radial grid, so we multiply again
        # by the grid for consistency.
        ϕ̃ = map(0:lmax) do l
            map(iχ_upf[l]) do i
                χln = upf.pswfc[i].chi
                χln = χln .* @view r[1:length(χln)]
                return χln
            end
        end
        ϕ̃ = OffsetVector(ϕ̃, 0:lmax)
    else
        ϕ̃ = nothing
    end

    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, D, ϕ̃, ρcore, ρval)
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

    # PSP8s store the projectors without any prefactor, so we multiply by the grid squared
    # for consistency.
    β = map(0:lmax) do l
        map(eachindex(psp8.projectors[l + 1])) do n
            βln = psp8.projectors[l + 1][n]
            βln = βln .* (@view r[1:length(βln)]).^2
            return βln  # _truncate(βln; atol=1e-12)
        end
    end
    β = OffsetVector(β, 0:lmax)

    D = OffsetVector(map(l -> diagm(psp8.ekb[l + 1]), 0:lmax), 0:lmax)
    ϕ̃ = nothing  # PSP8 doesn't support pseudo-atomic wavefunctions
    # PSP8s store the core charge density with a prefactor of 4π, which we remove for
    # consistency.
    ρcore = isnothing(psp8.rhoc) ? nothing : psp8.rhoc ./ 4π
    ρval = nothing  # PSP8 doesn't support pseudo-atomic valence charge density
    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, D, ϕ̃, ρcore, ρval)
end

is_norm_conserving(::NormConservingPsP)::Bool = true
is_ultrasoft(::NormConservingPsP)::Bool = false
is_paw(::NormConservingPsP)::Bool = false
