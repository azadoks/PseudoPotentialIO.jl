function get_upf_version(io::IO)::Int
    pos = position(io)
    seek(io, 0)
    line = readline(io)
    seek(io, pos)
    if occursin("<PP_INFO>", line)
        # Old UPF files start with the `<PP_INFO>` section
        return 1
    elseif occursin("UPF version=\"2.0.1\"", line)
        # New UPF files with schema are in XML and start with a version tag
        return 2
    else
        error("Unknown UPF version")
    end
end

function get_upf_version(path::AbstractString)::Int
    open(path, "r") do io
        return get_upf_version(io)
    end
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_HEADER>`
"""
struct UpfHeader
    "Generation code"
    generated::Union{Nothing,String}
    author::Union{Nothing,String}
    "Generation date"
    date::Union{Nothing,String}
    comment::Union{Nothing,String}
    "A valid chemical symbol: `{H, He, Li, ..., Og}`"
    element::String
    "A valid type of pseudopotential: `{NC, SL, 1/r, US, PAW, USPP}`"
    pseudo_type::String
    "A valid relativistic treatment: `{scalar, full, relativistic}`"
    relativistic::Union{Nothing,String}
    is_ultrasoft::Bool
    is_paw::Bool
    "True of the pseudopotential is just a bare Coulomb potential (all-electron)"
    is_coulomb::Bool
    "True if fully-relativistic with spin-orbit terms"
    has_so::Bool
    "True if all-electron wavefunctions present"
    has_wfc::Bool
    "True if data for GIPAW reconstruction is present"
    has_gipaw::Bool
    "True if data for GIPAW reconstruction is present"
    paw_as_gipaw::Union{Nothing,Bool}
    "True if non-linear core correction is included"
    core_correction::Bool
    "QuantumEspresso exchange-correlation identifier"
    functional::String
    "Pseudo-atomic charge"
    z_valence::Float64
    "Total pseudo-valence energy of the pseudopotential"
    total_psenergy::Union{Nothing,Float64}
    "Suggested plane wave cutoff for expansion of Kohn-Sham orbitals"
    wfc_cutoff::Union{Nothing,Float64}
    "Suggested plane wave cutoff for expansion of charge density"
    rho_cutoff::Union{Nothing,Float64}
    "Maximum angular momentum channel in the pseudopotential"
    l_max::Int
    "Maximum angular momentum channel in the atomic charge density (PAW only)"
    l_max_rho::Union{Nothing,Int}
    "Angular momentum chosen to be the local potential (-1 if none)"
    l_local::Union{Nothing,Int}
    "Number of points in the radial grid"
    mesh_size::Int
    "Number of pseudo-atomic wavefunctions"
    number_of_wfc::Int
    "Number of Kleinman-Bylander nonlocal projectors"
    number_of_proj::Int
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_MESH>`
"""
struct UpfMesh
    "Radial mesh"
    r::Vector{Float64}
    "Integration factors for integrating quantities on the radial mesh"
    rab::Vector{Float64}
    "Number of points in the radial mesh"
    mesh::Union{Nothing,Int}
    "Maximum value of the radial mesh"
    rmax::Union{Nothing,Float64}
    # Mesh generation parameters
    dx::Union{Nothing,Float64}
    xmin::Union{Nothing,Float64}
    zmesh::Union{Nothing,Float64}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJ.[i].[j]>`
"""
struct UpfQij
    qij::Vector{Float64}
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJL.[i].[j].[l]>`
"""
struct UpfQijl
    qijl::Vector{Float64}
    angular_momentum::Int
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_NONLOCAL/PP_AUGMENTATION>`
"""
struct UpfAugmentation
    """Integrals of the augmentation functions 4π ∫ Qij(r) r^2 dr.
    NB: `q = 0` does _not_ guarantee that the corresponding augmentation function is
    zero."""
    q::Matrix{Float64}
    """(PAW) Electronic multipoles of the corresponding augmentation channel. If the
    absolute value of a multipole is less than `augmentation_epsilon`, the corresponding
    augmentation function should be considered zero"""
    multipoles::Union{Nothing,Vector{Float64}}
    """Coefficients used to perform a Taylor expansion of the augmentation functions at
    small radii (NB: compulsory if `nqf > 0`, ignored otherwise)"""
    qfcoeff::Union{Nothing,Vector{Vector{Float64}}}
    rinner::Union{Nothing,Vector{Float64}}
    "If `q_with_l` is false, the augmentation functions for `i,j in 1:number_of_proj`"
    qijs::Union{Nothing,Vector{UpfQij}}
    """If `q_with_l` is true, the angular-momentum dependent augmentation functions for
    `i,j in 1:number_of_proj` and `l in 0:l_max_rho`"""
    qijls::Union{Nothing,Vector{UpfQijl}}
    "True if augmentation charge functions are decomposed into angular momentum components"
    q_with_l::Bool
    """Number of expansion coefficients for analytical expansion of the augmentation
    charge at small radius."""
    nqf::Int
    "Number of angular momenta terms in the augmentation charge, unused if `nqf = 0`"
    nqlc::Union{Nothing,Float64}
    """(UNUSED) (PAW) Shape of the augmentation function: `{PSQ, GAUSS, BESSEL}`, could
    be used for analyical overlap of PAW augmentation charge"""
    shape::Union{Nothing,String}
    "(PAW) Radial grid index beyond which augmentation charge is zero"
    iraug::Union{Nothing,Int}
    "(PAW) Radial distance beyond which augmentation charge is zero"
    raug::Union{Nothing,Float64}
    "(PAW): Maximum angular momentum appearing in augmentation charge"
    l_max_aug::Union{Nothing,Float64}
    """(PAW): Augmentation functions whose norms are less than `augmentation_epsilon` are
    considered zero"""
    augmentation_epsilon::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r_index::Union{Nothing,Float64}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_NONLOCAL/PP_BETA.[i]>`
"""
struct UpfBeta
    "Kleinman-Bylander nonlocal projector multiplied by the radial mesh, on the radial mesh"
    beta::Vector{Float64}
    index::Union{Nothing,Int}
    angular_momentum::Int
    cutoff_radius_index::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    norm_conserving_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    label::Union{Nothing,String}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_NONLOCAL>`
"""
struct UpfNonlocal
    """Kleinman-Bylander nonlocal projectors multiplied by the radial mesh,
    on the radial mesh"""
    betas::Vector{UpfBeta}
    "Kleinman-Bylander energies"
    dij::Matrix{Float64}
    "Agumentation data for ultrasoft and PAW pseudopotentials"
    augmentation::Union{Nothing,UpfAugmentation}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_PSWFC/PP_CHI>`
"""
struct UpfChi
    "Pseudo-atomic valence wavefunction on the radial mesh"
    chi::Vector{Float64}
    "Angular momentum"
    l::Int
    occupation::Float64
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    "Principal quantum number"
    n::Union{Nothing,Int}
    pseudo_energy::Union{Nothing,Float64}
    cutoff_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_SPIN_ORB/PP_RELWFC.[i]>`
"""
struct UpfRelWfc
    jchi::Float64
    index::Union{Nothing,Int}
    els::Union{Nothing,String}
    nn::Union{Nothing,Int}
    lchi::Union{Nothing,Int}
    oc::Union{Nothing,Float64}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_SPIN_ORB/PP_RELBETA.[i]>`
"""
struct UpfRelBeta
    index::Union{Nothing,Int}
    jjj::Float64
    lll::Union{Nothing,Int}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_SPIN_ORB>`
"""
struct UpfSpinOrb
    relwfcs::Vector{UpfRelWfc}
    relbetas::Vector{UpfRelBeta}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<//PP_(AE|PS)WFC.[i]>`
"""
struct UpfWfc
    wfc::Vector{Float64}
    index::Int
    l::Int
    label::Union{Nothing,String}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_FULL_WFC>`
"""
struct UpfFullWfc
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_PAW>`
"""
struct UpfPaw
    paw_data_format::Union{Nothing,Int}
    core_energy::Union{Nothing,Float64}
    occupations::Vector{Float64}
    ae_nlcc::Vector{Float64}
    ae_vloc::Vector{Float64}
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_GIPAW/PP_GIPAW_CORE_ORBITALS/PP_GIPAW_CORE_ORBITAL.[i]>`
"""
struct UpfGipawCoreOrbital
    index::Int
    label::Union{Nothing,String}
    "Principal quantum number"
    n::Int
    "Angular momentum"
    l::Int
    core_orbital::Vector{Float64}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF `<PP_GIPAW>`
"""
struct UpfGipaw
    gipaw_data_format::Int
    core_orbitals::Vector{UpfGipawCoreOrbital}
end

"""
$TYPEDEF

$(TYPEDFIELDS)

UPF pseudopotential
"""
struct UpfPsP <: PseudoPotentialIO.AbstractPsP
    "UPF format version"
    version::String
    "Optional general information about the pseudopotential, often generation input"
    info::Union{Nothing,String}
    "Various pseudopotential metadata"
    header::UpfHeader
    "Radial mesh, mesh integration factors, and other mesh information"
    mesh::UpfMesh
    "Pseudized core charge on the radial grid, (ignored if `core_correction` is false)"
    nlcc::Union{Nothing,Vector{Float64}}  # Σ_{i} 4π r_{i}^2 nlcc_{i}
    "Local part of the pseudopotential on the radial grid (ignored if `is_coulomb`)"
    local_::Union{Nothing,Vector{Float64}}
    "Nonlocal part of the pseudopotential"
    nonlocal::UpfNonlocal
    "Pseudo-atomic valence wavefunctions"
    pswfc::Union{Nothing,Vector{UpfChi}}
    "All-electron wavefunctions"
    full_wfc::Union{Nothing,UpfFullWfc}
    "Pseudo-atomic valence charge density on the radial grid"
    rhoatom::Vector{Float64}
    "Spin-orbit coupling data, (ignored if `has_so` is false)"
    spin_orb::Union{Nothing,UpfSpinOrb}
    "PAW data, (ignored if `is_paw` is false)"
    paw::Union{Nothing,UpfPaw}
    "GIPAW data"
    gipaw::Union{Nothing,UpfGipaw}
end

function UpfPsP(path::AbstractString)
    open(path, "r") do io
        return UpfPsP(io)
    end
end

function UpfPsP(io::IO)
    version = get_upf_version(io)
    if version == 1
        return upf1_parse_psp(io)
    end
    if version == 2
        text = read(io, String)
        # Remove end-of-file junk (input data, etc.)
        text = string(split(text, "</UPF>")[1], "</UPF>")
        # Clean any errant `&` characters
        text = replace(text, "&" => "")
        return upf2_parse_psp(parsexml(text))
    end
end

function element(psp::UpfPsP)::PeriodicTable.Element
    return PeriodicTable.elements[Symbol(psp.header.element)]
end
l_max(psp::UpfPsP)::Int = psp.header.l_max
n_proj_radial(psp::UpfPsP, l::Integer)::Int = count(beta -> beta.angular_momentum == l,
                                           psp.nonlocal.betas)
n_pseudo_wfc(psp::UpfPsP)::Int = psp.header.number_of_wfc
z_valence(psp::UpfPsP)::Float64 = psp.header.z_valence
is_paw(psp::UpfPsP)::Bool = psp.header.is_paw
is_ultrasoft(psp::UpfPsP)::Bool = psp.header.is_ultrasoft
is_norm_conserving(psp::UpfPsP)::Bool = psp.header.pseudo_type == "NC"
is_coulomb(psp::UpfPsP)::Bool = psp.header.is_coulomb
has_spin_orbit(psp::UpfPsP)::Bool = psp.header.has_so
has_nlcc(psp::UpfPsP)::Bool = psp.header.core_correction
relativistic_treatment(psp::UpfPsP)::Symbol = Symbol(psp.header.relativistic)
format(psp::UpfPsP)::String = "UPF v$(psp.version)"

function get_projector_radial(psp::UpfPsP, l, n)::UpfBeta
    projectors_l = filter(beta -> beta.angular_momentum == l, psp.nonlocal.betas)
    return projectors_l[n]
end

function get_pseudo_orbital_radial(psp::UpfPsP, l)::UpfChi
    orbitals = filter(orbital -> orbital.l == l, psp.pswfc)
    @assert length(orbital) == 1
    return first(orbitals)
end

function e_kb(psp::UpfPsP, li, ni, lj, nj)::Float64
    projector_li_ni = get_projector_radial(psp, li, ni)
    index_li_ni = projector_li_ni.index
    projector_lj_nj = get_projector_radial(psp, lj, nj)
    index_lj_nj = projector_lj_nj.index
    return psp.nonlocal.dij[index_li_ni, index_lj_nj]
end

function e_kb(psp::UpfPsP, l, n)::Float64
    projector_l_n = get_projector_radial(psp, l, n)
    index_l_n = projector_l_n.index
    return psp.nonlocal.dij[index_l_n, index_l_n]
end

function v_local_real(psp::UpfPsP, r::T)::T where {T<:Real}
    interpolator = linear_interpolation((psp.mesh.r,), psp.local_)
    return interpolator(r)
end

function v_local_fourier(psp::UpfPsP, q::T)::T where {T<:Real}
    v_corr_fourier = v_local_correction_fourier(psp, q)
    @. integrand = psp.mesh.r^2 * sphericalbesselj_fast(0, q * psp.mesh.r) *
                   (psp.local_ - v_local_correction_real(psp, psp.mesh.r))
    return 4T(π) * simpson(integrand, psp.rab) + v_corr_fourier
end

#TODO this will create NaNs for linear meshes that start at r=0
function projector_radial_real(psp::UpfPsP, l, n, r::T)::T where {T<:Real}
    projector = get_projector_radial(psp, l, n)
    interpolator = linear_interpolation((psp.mesh.r,), projector.beta)
    return interpolator(r) / r
end

function projector_radial_fourier(psp::UpfPsP, l, n, q::T)::T where {T<:Real}
    projector = get_projector_radial(psp, l, n)
    @. integrand = psp.mesh.r * sphericalbesselj_fast(l, q * psp.mesh.r) * projector
    return 4T(π) * trapezoid(integrand, psp.mesh.rab)
end

function pseudo_energy_correction(psp::UpfPsP)::Float64
    v_local_corrected = psp.local_ - v_local_correction_real.(psp, psp.mesh.r)
    return 4π * trapezoid(v_local_corrected, psp.mesh.rab)
end

function core_charge_density_real(psp::UpfPsP, r::T)::T where {T<:Real}
    interpolator = linear_interpolation((psp.mesh.r,), psp.nlcc)
    return interpolator(r)
end

function core_charge_density_fourier(psp::UpfPsP, q::T)::T where {T<:Real}
    @. integrand = psp.mesh.r^2 * sphericalbesselj_fast(0, q * psp.mesh.r) * psp.nlcc
    return 4T(π) * trapezoid(integrand, psp.mesh.rab)
end

#TODO this will create NaNs for linear meshes that start at r=0
function valence_charge_density_real(psp::UpfPsP, r::T)::T where {T<:Real}
    interpolator = linear_interpolation((psp.rhoatom,), psp.rhoatom)
    return interpolator(r) / (4T(π) * r^2)
end

function valence_charge_density_fourier(psp::UpfPsP, q::T)::T where {T<:Real}
    @. integrand = sphericalbesselj_fast(0, q * psp.mesh.r) * psp.rhoatom
    return trapezoid(integrand, psp.mesh.rab)
end

function pseudo_orbital_radial_real(psp::UpfPsP, l, r::T)::T where {T<:Real}
    orbital = get_pseudo_orbital_radial(psp, l)
    interpolator = linear_interpolation((psp.rhoatom,), orbital.chi)
    return interpolator(r)
end

function pseudo_orbital_radial_fourier(psp::UpfPsP, l, q::T)::T where {T<:Real}
    orbital = get_pseudo_orbital_radial(psp, l)
    @. integrand = psp.mesh.r^2 * sphericalbesselj_fast(0, q * psp.mesh.r) * orbital
    return 4T(π) * trapezoid(integrand, psp.mesh.rab)
end
