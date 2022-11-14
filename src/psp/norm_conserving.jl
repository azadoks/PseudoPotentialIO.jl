"""
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
    ρcore = isnothing(upf.nlcc) ? nothing : _truncate(upf.nlcc; atol=1e-8)

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

    # Find the cutoff radius index for each projector
    β_ircut = map(0:lmax) do l
        map(iβ_upf[l]) do i
            return length(upf.nonlocal.betas[i].beta)
        end
    end
    β_ircut = OffsetVector(β_ircut, 0:lmax)

    # UPFs store the projectors multiplied by the radial grid. For compatability with the
    # PseudoPotentialIO interface, we need the pure projectors. To get them, we divide
    # rβ by r where r > 0 and use that data to extrapolate to r = 0 using a quadratic
    # polynomial.
    β = map(0:lmax) do l
        map(iβ_upf[l]) do i
            return _standardize_wavefunction_like(upf.nonlocal.betas[i].beta, r)
        end
    end
    β = OffsetVector(β, 0:lmax) ./ 2  # Ry -> Ha

    # UPFs store the pseudo-atomic valence charge density with a prefactor of 4π r^2.
    # We do the extrapolation like for the projectors.
    ρval = _standardize_ρval(upf.rhoatom, r)
    ρval = _truncate(ρval; atol=1e-8)

    # The pseudo-atomic wavefunctions need the same extrapolation treatment as the
    # projectors.
    if !isnothing(upf.pswfc)
        # Collect the indices in upf.nonlocal.betas for projectors at each angular momentum
        iχ_upf = map(0:lmax) do l
            return filter(i -> upf.pswfc[i].l == l, eachindex(upf.pswfc))
        end
        iχ_upf = OffsetVector(iχ_upf, 0:lmax)

        ϕ̃ = map(0:lmax) do l
            map(iχ_upf[l]) do i
                χ_std = _standardize_wavefunction_like(upf.pswfc[i].chi, r)
                χ_std = _truncate(χ_std; atol=1e-4)
            end
        end
        ϕ̃ = OffsetVector(ϕ̃, 0:lmax)

        ϕ̃_ircut = OffsetVector(map(l -> length.(ϕ̃[l]), 0:lmax), 0:lmax)
    else
        ϕ̃_ircut = nothing
        ϕ̃ = nothing
    end

    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, β_ircut, D, ϕ̃,
                                      ϕ̃_ircut, ρcore, ρval)
end

function _standardize_ρval(ρval_upf::Vector{Float64}, r::Vector{Float64})::Vector{Float64}
    ρval = @. ρval_upf / (4Float64(π) * r^2)
    if iszero(r[1])
        ρval[1] = _extrapolate_to_zero(ρval, r)
    end
    return ρval
end

function _standardize_wavefunction_like(ϕ_upf::Vector{Float64},
                                        r::Vector{Float64})::Vector{Float64}
    irmax = length(ϕ_upf)
    ϕ = ϕ_upf[1:irmax] ./ r[1:irmax]
    if iszero(r[1])
        ϕ[1] = _extrapolate_to_zero(ϕ, r)
    end
    return ϕ
end

function _extrapolate_to_zero(f::AbstractVector{Float64}, r::AbstractVector{Float64};
                              atol::Float64=1e-1)::Float64
    val = fit(r[2:6], f[2:6], 2)(r[1])
    # If f is approaching 0 at r = 0, force it to numerical 0
    return isapprox(val, 0; atol) ? 0.0 : val
end

function _truncate(f::AbstractVector{Float64}; atol=sqrt(eps(Float64)), length_min=2)
    icut = findfirst(i -> maximum(abs, @view f[i:end]; init=0.) < atol, eachindex(f))
    icut = isnothing(icut) ? firstindex(f) + length_min - 1 : icut
    icut = (icut - firstindex(f)) < length_min ? firstindex(f) + length_min - 1 : icut
    return f[firstindex(f):icut]
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

    # Cut off the projectors at the first index after which their mean value goes to
    # approximately zero.
    β = map(0:lmax) do l
        map(eachindex(psp8.projectors[l + 1])) do i
            proj = psp8.projectors[l + 1][i]
            proj = _truncate(proj; atol=1e-8)
        end
    end
    β = OffsetVector(β, 0:lmax)

    β_ircut = OffsetVector(map(l -> length.(β[l]), 0:lmax), 0:lmax)
    D = OffsetVector(map(l -> diagm(psp8.ekb[l + 1]), 0:lmax), 0:lmax)
    ϕ̃ = nothing  # PSP8 doesn't support pseudo-atomic wavefunctions
    ϕ̃_ircut = nothing
    ρcore = psp8.rhoc ./ (4π)
    ρval = nothing  # PSP8 doesn't support pseudo-atomic valence charge density
    return NormConservingPsP{Float64}(Ztot, Zval, lmax, r, dr, Vloc, β, β_ircut, D, ϕ̃,
                                      ϕ̃_ircut, ρcore, ρval)
end

is_norm_conserving(::NormConservingPsP)::Bool = true
is_ultrasoft(::NormConservingPsP)::Bool = false
is_paw(::NormConservingPsP)::Bool = false
