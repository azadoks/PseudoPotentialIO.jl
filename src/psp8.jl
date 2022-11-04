function parse_fortran(::Type{T}, x::AbstractString) where {T<:Real}
    return parse(T, replace(lowercase(x), "d" => "e"))
end

"""
$TYPEDEF

$(TYPEDFIELDS)

PSP8 header block
"""
struct Psp8Header
    "Description"
    title::String
    "Atomic charge"
    zatom::Float64
    "Pseudo-atomic (valence) charge"
    zion::Float64
    "Generation date `ddmmyy`"
    pspd::Int
    "PSP format version"
    pspcod::Int
    "ABINIT code for exchange-correlation if positive, LibXC code if negative"
    pspxc::Int
    "Maximum angular momentum channel of the Kleinman-Bylander projectors"
    lmax::Int
    "Angular momentum channel of the local potential (`>lmax` if none)"
    lloc::Int
    "Number of points on the radial mesh"
    mmax::Int
    "Unused"
    r2well::Float64
    "Radius beyond which the model core charge (if present) is zero or negligible"
    rchrg::Float64
    "Pseudopotential has a model core charge if > 0"
    fchrg::Float64
    "Unused"
    qchrg::Float64
    "Number of projectors for each angular momentum channel"
    nproj::Vector{Int}
    "Signals presence of spin-orbit coupling if 2 or 3"
    extension_switch::Int
    "Number of spin-orbit projectors for each angular momentum (if present)"
    nprojso::Union{Nothing,Vector{Int}}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

PSP8 pseudopotential
"""
struct Psp8PsP <: PseudoPotentialIO.AbstractPsP
    "Various pseudopotential metadata"
    header::Psp8Header
    "Radial grid"
    rgrid::Vector{Float64}
    "Local part of the pseudopotential"
    v_local::Vector{Float64}
    "Radial part of the Kleinman-Bylander projectors for each angular momentum"
    projectors::Vector{Vector{Vector{Float64}}}
    "Kleinman-Bylander energies for each angular momentum"
    ekb::Vector{Vector{Float64}}
    "Radial part of the spin-orbit Kleinman-Bylander projectors for each angular momentum"
    projectors_so::Union{Nothing,Vector{Vector{Vector{Float64}}}}
    "Spin-orbit Kleinman-Bylander energies for each angular momentum"
    ekb_so::Union{Nothing,Vector{Vector{Float64}}}
    "Model core charge"
    rhoc::Union{Nothing,Vector{Float64}}
    "First derivative of the model core charge"
    d_rhoc_dr::Union{Nothing,Vector{Float64}}
    "Second derivative of the model core charge"
    d2_rhoc_dr2::Union{Nothing,Vector{Float64}}
    "Third derivative of the model core charge"
    d3_rhoc_dr3::Union{Nothing,Vector{Float64}}
    "Fourth derivative of the model core charge"
    d4_rhoc_dr4::Union{Nothing,Vector{Float64}}
end

function psp8_parse_header(io::IO)
    # line 1
    title = readline(io)
    # line 2
    s = split(readline(io))
    zatom = parse_fortran(Float64, s[1])
    zion = parse_fortran(Float64, s[2])
    pspd = parse(Int, s[3])
    # line 3
    s = split(readline(io))
    pspcod = parse(Int, s[1])
    pspxc = parse(Int, s[2])
    lmax = parse(Int, s[3])
    lloc = parse(Int, s[4])
    mmax = parse(Int, s[5])
    r2well = parse(Int, s[6])
    @assert pspcod == 8
    # line 4
    s = split(readline(io))
    rchrg = parse_fortran(Float64, s[1])
    fchrg = parse_fortran(Float64, s[2])
    qchrg = parse_fortran(Float64, s[3])
    # line 5
    s = split(readline(io))
    nproj = [parse(Int, s[l + 1]) for l in 0:lmax]
    if lloc <= lmax
        @assert nproj[lloc + 1] == 0
    end
    # line 6
    s = split(readline(io))
    extension_switch = parse(Int, s[1])
    # line 7
    if extension_switch in (2, 3)
        s = split(readline(io))
        nprojso = [0, [parse(Int, s[i]) for i in 1:lmax]...]
    else
        nprojso = nothing
    end

    return Psp8Header(title, zatom, zion, pspd, pspcod, pspxc, lmax, lloc, mmax, r2well,
                      rchrg, fchrg, qchrg, nproj, extension_switch, nprojso)
end

function psp8_parse_projector_block(io, nproj, mmax)
    header_line = split(readline(io))
    l = parse(Int, header_line[1])
    nproj_l = nproj[l + 1]
    ekb = [parse_fortran(Float64, header_line[i + 1]) for i in 1:nproj_l]

    rgrid = Vector{Float64}(undef, mmax)
    projectors = [Vector{Float64}(undef, mmax) for _ in 1:(nproj_l)]
    for i in 1:mmax
        s = split(readline(io))
        rgrid[i] = parse_fortran(Float64, s[2])
        for j in 1:nproj_l
            projectors[j][i] = parse_fortran(Float64, s[2 + j])
        end
    end

    return (; l, rgrid, projectors, ekb)
end

function psp8_parse_v_local(io, mmax)
    header_line = split(readline(io))
    l = parse(Int, header_line[1])

    rgrid = Vector{Float64}(undef, mmax)
    v_local = Vector{Float64}(undef, mmax)
    for i in 1:mmax
        s = split(readline(io))
        rgrid[i] = parse_fortran(Float64, s[2])
        v_local[i] = parse_fortran(Float64, s[3])
    end

    return (; l, rgrid, v_local)
end

function psp8_parse_main_blocks(io, mmax, nproj, lmax, lloc)
    projector_blocks = []
    v_local_block = ()

    if lmax < lloc
        n_blocks = lmax + 1
    else
        n_blocks = lmax
    end

    for _ in 0:n_blocks
        # Record the position at the start of the block so we can
        # read in the first line and go back
        pos = position(io)
        # Read the block header
        header_line = split(readline(io))
        # Go back to the start of the block
        seek(io, pos)
        # Parse the block's angular momentum
        block_l = parse(Int, header_line[1])
        if block_l == lloc
            v_local_block = psp8_parse_v_local(io, mmax)
            if lloc <= lmax
                push!(projector_blocks, (; l=block_l, rgrid=nothing, projectors=[], ekb=[]))
            end
        else
            block = psp8_parse_projector_block(io, nproj, mmax)
            push!(projector_blocks, block)
        end
    end

    projectors = [block.projectors for block in projector_blocks]
    ekb = [block.ekb for block in projector_blocks]

    return (; v_local_block.rgrid, v_local_block.v_local, projectors, ekb)
end

function psp8_parse_spin_orbit_blocks(io, mmax, nprojso, lmax)
    projector_blocks = []
    for _ in 1:lmax
        block = psp8_parse_projector_block(io, nprojso, mmax)
        push!(projector_blocks, block)
    end

    projectors = [block.projectors for block in projector_blocks]
    ekb = [block.ekb for block in projector_blocks]

    return (; projectors, ekb)
end

function psp8_parse_nlcc(io, mmax)
    radial_grid = Vector{Float64}(undef, mmax)
    rhoc = Vector{Float64}(undef, mmax)
    d_rhoc_dr = Vector{Float64}(undef, mmax)
    d2_rhoc_dr2 = Vector{Float64}(undef, mmax)
    d3_rhoc_dr3 = Vector{Float64}(undef, mmax)
    d4_rhoc_dr4 = Vector{Float64}(undef, mmax)
    for i in 1:mmax
        s = split(readline(io))
        radial_grid[i] = parse_fortran(Float64, s[2])
        # These include the 4π factor
        rhoc[i] = parse_fortran(Float64, s[3])
        d_rhoc_dr[i] = parse_fortran(Float64, s[4])
        d2_rhoc_dr2[i] = parse_fortran(Float64, s[5])
        d3_rhoc_dr3[i] = parse_fortran(Float64, s[6])
        d4_rhoc_dr4[i] = parse_fortran(Float64, s[7])
    end
    return (; rhoc, d_rhoc_dr, d2_rhoc_dr2, d3_rhoc_dr3, d4_rhoc_dr4)
end

function Psp8PsP(io::IO)
    header = psp8_parse_header(io)
    main_blocks = psp8_parse_main_blocks(io, header.mmax, header.nproj, header.lmax,
                                         header.lloc)
    if header.extension_switch in (2, 3)
        spin_orbit = psp8_parse_spin_orbit_blocks(io, header.mmax, header.nprojso,
                                                  header.lmax)
        # Add empty vectors for l=0 to maintain angular momentum - index correspondence
        spin_orbit = (projectors=[[], spin_orbit.projectors...],
                      ekb=[[], spin_orbit.ekb...])
    else
        spin_orbit = (projectors=nothing, ekb=nothing)
    end
    if header.fchrg > 0
        nlcc = psp8_parse_nlcc(io, header.mmax)
    else
        nlcc = (rhoc=nothing, d_rhoc_dr=nothing, d2_rhoc_dr2=nothing, d3_rhoc_dr3=nothing,
                d4_rhoc_dr4=nothing)
    end
    return Psp8PsP(header, main_blocks.rgrid, main_blocks.v_local, main_blocks.projectors,
                   main_blocks.ekb, spin_orbit.projectors, spin_orbit.ekb,
                   nlcc.rhoc, nlcc.d_rhoc_dr, nlcc.d2_rhoc_dr2, nlcc.d3_rhoc_dr3,
                   nlcc.d4_rhoc_dr4)
end

function Psp8PsP(path::AbstractString)
    open(path, "r") do io
        return Psp8PsP(io)
    end
end

function element(psp::Psp8PsP)::PeriodicTable.Element
    return PeriodicTable.elements[round(Int, psp.header.zatom)]
end
l_max(psp::Psp8PsP)::Int = psp.header.lmax
n_proj_radial(psp::Psp8PsP, l::Integer)::Int = psp.header.nproj[l + 1]
n_pseudo_wfc(::Psp8PsP)::Int = 0
z_valence(psp::Psp8PsP)::Float64 = psp.header.zion
is_paw(::Psp8PsP)::Bool = false
is_ultrasoft(::Psp8PsP)::Bool = false
is_norm_conserving(::Psp8PsP)::Bool = true
is_coulomb(::Psp8PsP)::Bool = false
has_spin_orbit(psp::Psp8PsP)::Bool = psp.header.extension_switch in (2, 3)
has_nlcc(psp::Psp8PsP)::Bool = psp.header.fchrg > 0
relativistic_treatment(psp::Psp8PsP)::Symbol = has_spin_orbit(psp) ? :full : :scalar
formalism(::Psp8PsP)::Symbol = :norm_conserving
format(::Psp8PsP) = "PSP v8"
dr(psp::Psp8PsP) = mean(diff(psp.rgrid))

function get_projector_radial(psp::Psp8PsP, l::Integer, n::Integer)::Vector{Float64}
    return psp.projectors[l + 1][n]
end

function e_kb(psp::Psp8PsP, l::Integer, n::Integer)::Vector{Float64}
    return psp.ekb[l][n]
end

function v_local_real(psp::Psp8PsP, r::T)::T where {T<:Real}
    interpolator = linear_interpolation((psp.rgrid,), psp.v_local)
    return interpolator(r)
end

function v_local_fourier(psp::Psp8PsP, q::T)::T where {T<:Real}
    v_corr_fourier = v_local_correction_fourier(psp, q)
    @. integrand = psp.rgrid^2 * sphericalbesselj_fast(0, q * psp.rgrid) *
                   (psp.v_local - v_local_correction_real(psp, psp.rgrid))
    return 4T(π) * simpson(integrand, dr(psp)) + v_corr_fourier
end

function projector_radial_real(psp::Psp8PsP, l::Integer, n::Integer,
                               r::T)::T where {T<:Real}
    projector = get_projector_radial(psp, l, n)
    interpolator = linear_interpolation((psp.rgrid,), projector)
    return interpolator(r)
end

function projector_radial_fourier(psp::Psp8PsP, l::Integer, n::Integer,
                                  q::T)::T where {T<:Real}
    projector = get_projector_radial(psp, l, n)
    @. integrand = psp.rgrid^2 * sphericalbesselj_fast(l, q * psp.rgrid) * projector
    return 4T(π) * trapezoid(integrand, dr(psp))
end

function pseudo_energy_correction(psp::Psp8PsP)::Float64
    v_local_corrected = psp.v_local - v_local_correction_real.(psp, psp.rgrid)
    return 4Float64(π) * trapezoid(v_local_corrected, dr(psp))
end

function core_charge_density_real(psp::Psp8PsP, r::T)::T where {T<:Real}
    interpolator = linear_interpolation((psp.rgrid,), psp.rhoc)
    return interpolator(r)
end

function core_charge_density_fourier(psp::Psp8PsP, q::T)::T where {T<:Real}
    @. integrand = psp.rgrid^2 * sphericalbesselj_fast(0, q * psp.rgrid) * psp.rhoc
    return trapezoid(integrand, dr(psp))
end
