@doc raw"""
Type representing a numeric norm-conserving pseudopotential.
"""
struct NormConservingPsP{T} <: NumericPsP{T}
    "Identifier"
    identifier::String
    "SHA1 Checksum"
    checksum::Union{Nothing,Vector{UInt8}}
    "Total charge."
    Zatom::Int
    "Valence charge."
    Zval::Int
    "Maximum angular momentum."
    lmax::Int
    "Local part of the potential on the radial mesh (without an r² prefactor)."
    Vloc::NumericNumericLocalPotential{T}
    "Nonlocal projectors β[l][n] on the radial mesh (with an r² prefactor)."
    β::OffsetVector{Vector{NumericProjector{T}},Vector{Vector{NumericProjector{T}}}}
    "Projector coupling coefficients D[l][n,m]."
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions χ[l][n] on the radial mesh (with an r² prefactor)."
    χ::OffsetVector{Vector{NumericState{T}},Vector{Vector{NumericState{T}}}}
    "Model core charge density on the radial mesh (with an r² prefactor)."
    ρcore::Union{Nothing,NumericDensity{T}}
    "Valence charge density on the radial mesh (with an r² prefactor)."
    ρval::Union{Nothing,NumericDensity{T}}
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
    # The quantity that's used in calculations is ⟨β|D|β⟩, so the units don't practically
    # matter. However, HGH pseudos in UPF format use the first units, so we assume them
    # to facilitate comparison of the intermediate quantities with analytical HGH.

    Zatom = PeriodicTable.elements[Symbol(upf.header.element)].number
    Zval = upf.header.z_valence
    lmax = upf.header.l_max
    r = upf.mesh.r
    dr = upf.mesh.rab

    Vloc_data = upf.local_ ./ 2  # Ry -> Ha
    Vloc_mesh = ArbitraryMesh(r[eachindex(Vloc_data)], dr[eachindex(Vloc_data)])
    Vloc = NumericNumericLocalPotential(Vloc_mesh, Vloc_data)

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

    # UPFs store the projectors multiplied by the radial grid, so we multiply again by the
    # grid for consistency.
    # If the cutoff radius index is stored, we use it to truncate the projector for
    # numerical stability.
    β = map(0:lmax) do l
        map(iβ_upf[l]) do n
            βln_data = upf.nonlocal.betas[n].beta ./ 2  # Ry -> Ha

            ir_cut = upf.nonlocal.betas[n].cutoff_radius_index
            ir_cut = isnothing(ir_cut) ? length(βln) : ir_cut

            βln_data = @view(βln[1:ircut]) .* @view(r[1:ir_cut])  # rβln -> r²βln
            βln_mesh = ArbitraryMesh(r[eachindex(βln_data)], dr[eachindex(βln_data)])
            βln_jₗ = fast_sphericalbesselj(l)

            return NumericProjector(n, l, 0.0, βln_mesh, βln_data, βln_jₗ)
        end
    end
    β = OffsetVector(β, 0:lmax)

    # Extract the blocks from `upf.nonlocal.dij` corresponding to each angular momentum
    D = map(1:(length(cumul_nβ) - 1)) do i
        return collect(upf.nonlocal.dij[(cumul_nβ[i] + 1):cumul_nβ[i + 1],
                                        (cumul_nβ[i] + 1):cumul_nβ[i + 1]])
    end
    D = OffsetVector(D, 0:lmax) .* 2  # 1/Ry -> 1/Ha

    # UPFs store the pseudo-atomic valence charge density with a prefactor of 4πr².
    # For consistency, we remove the 4π prefactor.
    if isnothing(upf.rhoatom)
        ρval = nothing
    else
        ρval_data = upf.rhoatom ./ 4π
        ρval_mesh = ArbitraryMesh(r[eachindex(ρval_data)], dr[eachindex(ρval_data)])
        ρval = NumericDensity(ρval_mesh, ρval_data)
    end

    # UPFs store the core charge density as a true charge (without 4πr² as a prefactor),
    # so we multiply by r² for consistency
    if isnothing(upf.nlcc)
        ρcore = nothing
    else
        ρcore_mesh = ArbitraryMesh(r[eachindex(upf.nlcc)], dr[eachindex(upf.nlcc)])
        ρcore_data = upf.nlcc .* (r[eachindex(upf.nlcc)]) .^ 2
        ρcore = NumericDensity(ρcore_mesh, ρcore_data)
    end

    if !isnothing(upf.pswfc)
        # Collect the indices in upf.pswfc for projectors at each angular momentum
        iχ_upf = map(0:lmax) do l
            return filter(i -> upf.pswfc[i].l == l, eachindex(upf.pswfc))
        end
        iχ_upf = OffsetVector(iχ_upf, 0:lmax)

        # UPFs store the wavefunctions multiplied by the radial grid, so we multiply again
        # by r (to get r²χ) for consistency.
        χ = map(0:lmax) do l
            map(iχ_upf[l]) do i
                χln_data = upf.pswfc[i].chi
                χln_data = χln .* r[eachindex(χln_data)]
                χln_mesh = ArbitraryMesh(r[eachindex(χln_data)], dr[eachindex(χln_data)])
                χln_jₗ = fast_sphericalbesselj(l)
                return NumericState(n, l, 0.0, χln_mesh, χln_data, χln_jₗ)
            end
        end
        χ = OffsetVector(χ, 0:lmax)
    else
        χ = OffsetVector([NumericState[] for _ in 0:lmax])
    end

    return NormConservingPsP{Float64}(upf.identifier, upf.checksum, Zatom, Zval, lmax, r,
                                      Vloc, β, D, χ, ρcore, ρval)
end

function NormConservingPsP(psp8::Psp8File)
    if psp8.header.extension_switch in (2, 3)
        error("Fully relativistic pseudos are not supported")
    end

    Zatom = psp8.header.zatom
    Zval = psp8.header.zion
    lmax = psp8.header.lmax

    r = ArbitraryMesh(psp8.rgrid)

    Vloc = NumericNumericLocalPotential(r, psp8.v_local)

    # PSP8s store the projectors without any prefactor, so we multiply by r² for consitency
    β = map(0:lmax) do l
        map(eachindex(psp8.projectors[l + 1])) do n
            βln_data = psp8.projectors[l + 1][n]
            βln_data = βln .* r .^ 2
            βln_jₗ = fast_sphericalbesselj(l)
            βln = NumericProjector(n, l, 0.0, r, βln_data, βln_jₗ)
            return βln
        end
    end
    β = OffsetVector(β, 0:lmax)

    D = OffsetVector(map(l -> diagm(psp8.ekb[l + 1]), 0:lmax), 0:lmax)
    χ = OffsetVector([NumericState[] for _ in 0:lmax])  # PSP8 doesn't support chi-functions
    # PSP8s store the core charge density with a prefactor of 4π, so we remove it and
    # multiply by r² for consistency.
    ρcore = isnothing(psp8.rhoc) ? nothing : NumericDensity(r, psp8.rhoc .* r .^ 2 ./ 4π)
    ρval = isnothing(psp8.rhov) ? nothing : NumericDensity(r, psp8.rhov .* r .^ 2 ./ 4π)
    return NormConservingPsP{Float64}(psp8.identifier, psp8.checksum, Zatom, Zval, lmax, r,
                                      Vloc, β, D, χ, ρcore, ρval)
end

is_norm_conserving(::NormConservingPsP)::Bool = true
is_ultrasoft(::NormConservingPsP)::Bool = false
is_paw(::NormConservingPsP)::Bool = false
