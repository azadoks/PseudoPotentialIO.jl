"""
UPF `<PP_HEADER>`
"""
struct UpfHeader <: PsPFile
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
    "True if pseudo-atomic wavefunctions present"
    has_wfc::Bool
    "True if data for GIPAW reconstruction is present"
    has_gipaw::Bool
    "True if data for GIPAW reconstruction is present"
    paw_as_gipaw::Union{Nothing,Bool}
    "True if non-linear core correction is included"
    core_correction::Bool
    """
    True if additional data (model and atom kinetic energy densities)
    for kinetic energy-dependent meta-GGAs is present.
    """
    with_metagga_info::Bool
    "QuantumEspresso exchange-correlation identifiers"
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
    "Number of chi-functions"
    number_of_wfc::Int
    "Number of Kleinman-Bylander nonlocal projectors"
    number_of_proj::Int
    UpfHeader(fields...) = validated(new(fields...))
end

function Base.show(io::IO, ::MIME"text/plain", header::UpfHeader)
    println(io, typeof(header))
    if !isnothing(header.generated)
        @printf "%032s: %s\n" "generated" header.generated
    end
    if !isnothing(header.author)
        @printf "%032s: %s\n" "author" header.author
    end
    if !isnothing(header.date)
        @printf "%032s: %s\n" "date" header.date
    end
    @printf "%032s: %s\n" "element" header.element
    @printf "%032s: %s\n" "pseudopotential type" header.pseudo_type
    if !isnothing(header.relativistic)
        @printf "%032s: %s\n" "relativistic treatment" header.relativistic
    end
    @printf "%032s: %s\n" "is ultrasoft" header.is_ultrasoft
    @printf "%032s: %s\n" "is PAW" header.is_paw
    @printf "%032s: %s\n" "is bare Coulomb" header.is_coulomb
    @printf "%032s: %s\n" "has spin-orbit" header.has_so
    @printf "%032s: %s\n" "has pseudo-atomic wavefunctions" header.has_wfc
    @printf "%032s: %s\n" "has GIPAW data" header.has_gipaw
    if !isnothing(header.paw_as_gipaw)
        @printf "%032s: %s\n" "PAW as GIPAW" header.paw_as_gipaw
    end
    @printf "%032s: %s\n" "has non-linear core correction" header.core_correction
    @printf "%032s: %s\n" "has meta-GGA info" header.with_metagga_info
    @printf "%032s: %s\n" "functional (QuantumESPRESSO)" header.functional
    @printf "%032s: %.6f\n" "valence charge" header.z_valence
    if !isnothing(header.total_psenergy)
        @printf "%032s: %.6f\n" "total pseudo energy" header.total_psenergy
    end
    if !isnothing(header.wfc_cutoff)
        @printf "%032s: %.6f Ry\n" "suggested wavefunction cutoff" header.wfc_cutoff
    end
    if !isnothing(header.rho_cutoff)
        @printf "%032s: %.6f Ry\n" "suggested charge density cutoff" header.rho_cutoff
    end
    @printf "%032s: %d\n" "max angular momentum" header.l_max
    if !isnothing(header.l_max_rho)
        @printf "%032s: %d\n" "max angular momentum in charge density (PAW)" header.l_max_rho
    end
    @printf "%032s: %d\n" "local potential angular momentum" header.l_local
    @printf "%032s: %d\n" "radial mesh size" header.mesh_size
    @printf "%032s: %d\n" "number of pseudo wavefunctions" header.number_of_wfc
    @printf "%032s: %d\n" "number of KB projectors" header.number_of_proj
end

"""
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
    UpfMesh(fields...) = validated(new(fields...))
end

"""
UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJ.[i].[j]>`
"""
struct UpfQij
    qij::Vector{Float64}
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
    UpfQij(fields...) = validated(new(fields...))
end

"""
UPF `<PP_NONLOCAL/PP_AUGMENTATION/PP_QIJL.[i].[j].[l]>`
"""
struct UpfQijl
    qijl::Vector{Float64}
    angular_momentum::Int
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    composite_index::Union{Nothing,Int}
    is_null::Union{Nothing,Bool}
    UpfQijl(fields...) = validated(new(fields...))
end

struct UpfQfcoef
    qfcoef::Vector{Float64}
    angular_momentum::Union{Nothing,Int}
    first_index::Union{Nothing,Int}
    second_index::Union{Nothing,Int}
    UpfQfcoef(fields...) = validated(new(fields...))
end

"""
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
    qfcoefs::Union{Nothing,Vector{UpfQfcoef}}
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
    nqlc::Union{Nothing,Int}
    """(UNUSED) (PAW) Shape of the augmentation function: `{PSQ, GAUSS, BESSEL}`, could
    be used for analyical overlap of PAW augmentation charge"""
    shape::Union{Nothing,String}
    "(PAW) Radial grid index beyond which augmentation charge is zero"
    iraug::Union{Nothing,Int}
    "(PAW) Radial distance beyond which augmentation charge is zero"
    raug::Union{Nothing,Float64}
    "(PAW): Maximum angular momentum appearing in augmentation charge"
    l_max_aug::Union{Nothing,Int}
    """(PAW): Augmentation functions whose norms are less than `augmentation_epsilon` are
    considered zero"""
    augmentation_epsilon::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r::Union{Nothing,Float64}
    "(DEPRECATED?)"
    cutoff_r_index::Union{Nothing,Int}
    UpfAugmentation(fields...) = validated(new(fields...))
end

"""
UPF `<PP_NONLOCAL/PP_BETA.[i]>`
"""
struct UpfBeta
    "Kleinman-Bylander nonlocal projector multiplied by the radial mesh, on the radial mesh"
    beta::Vector{Float64}  # misnomer: actually r*β(r)!
    index::Union{Nothing,Int}
    angular_momentum::Int
    cutoff_radius_index::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    norm_conserving_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    label::Union{Nothing,String}
    UpfBeta(fields...) = validated(new(fields...))
end

"""
UPF `<PP_NONLOCAL>`
"""
struct UpfNonlocal
    """Kleinman-Bylander nonlocal projectors"""
    betas::Vector{UpfBeta}
    "Kleinman-Bylander energies"
    dij::Matrix{Float64}
    "Agumentation data for ultrasoft and PAW pseudopotentials"
    augmentation::Union{Nothing,UpfAugmentation}
    UpfNonlocal(fields...) = validated(new(fields...))
end

"""
UPF `<PP_PSWFC/PP_CHI>`
"""
struct UpfChi
    "Pseudo-atomic valence wavefunction on the radial mesh"
    chi::Vector{Float64}  # misnomer: actually r*χ(r)!
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
    UpfChi(fields...) = validated(new(fields...))
end

"""
UPF `<PP_SPIN_ORB/PP_RELWFC.[i]>`
"""
struct UpfRelWfc
    jchi::Float64
    index::Union{Nothing,Int}
    els::Union{Nothing,String}
    nn::Union{Nothing,Int}
    lchi::Union{Nothing,Int}
    oc::Union{Nothing,Float64}
    UpfRelWfc(fields...) = validated(new(fields...))
end

"""
UPF `<PP_SPIN_ORB/PP_RELBETA.[i]>`
"""
struct UpfRelBeta
    index::Union{Nothing,Int}
    jjj::Float64
    lll::Union{Nothing,Int}
    UpfRelBeta(fields...) = validated(new(fields...))
end

"""
UPF `<PP_SPIN_ORB>`
"""
struct UpfSpinOrb
    relwfcs::Vector{UpfRelWfc}
    relbetas::Vector{UpfRelBeta}
    UpfSpinOrb(fields...) = validated(new(fields...))
end

"""
UPF `<//PP_(AE|PS)WFC.[i]>`
"""
struct UpfWfc
    wfc::Vector{Float64}
    index::Int
    l::Int
    label::Union{Nothing,String}
    UpfWfc(fields...) = validated(new(fields...))
end

"""
UPF `<PP_FULL_WFC>`
"""
struct UpfFullWfc
    aewfcs::Vector{UpfWfc}
    pswfcs::Vector{UpfWfc}
    UpfFullWfc(fields...) = validated(new(fields...))
end

"""
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
    UpfPaw(fields...) = validated(new(fields...))
end

"""
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
    UpfGipawCoreOrbital(fields...) = validated(new(fields...))
end

"""
UPF `<PP_GIPAW>`
"""
struct UpfGipaw
    gipaw_data_format::Int
    core_orbitals::Vector{UpfGipawCoreOrbital}
    UpfGipaw(fields...) = validated(new(fields...))
end

"""
Universal Pseudopotential Format file contents.
"""
struct UpfFile <: PsPFile
    "Identifier"
    identifier::String
    "UPF format version"
    version::String
    "Optional general information about the pseudopotential, often generation input"
    info::Union{Nothing,String}
    "Various pseudopotential metadata"
    header::UpfHeader
    "Radial mesh, mesh integration factors, and other mesh information"
    mesh::UpfMesh
    """
    Pseudized core charge on the radial grid (ignored if `core_correction` is false).
    Note that there is no 4π nor r² prefactor here.
    """
    nlcc::Union{Nothing,Vector{Float64}}
    """
    Model core kinetic energy density on the radial grid.
    Note that there is no 4π nor r² prefactor here.
    """
    taumod::Union{Nothing,Vector{Float64}}
    "Local part of the pseudopotential on the radial grid (ignored if `is_coulomb`)"
    local_::Union{Nothing,Vector{Float64}}
    "Nonlocal part of the pseudopotential"
    nonlocal::UpfNonlocal
    "Pseudo-atomic valence wavefunctions"
    pswfc::Vector{UpfChi}
    "All-electron wavefunctions"
    full_wfc::Union{Nothing,UpfFullWfc}
    """
    Pseudo-atomic valence charge density on the radial grid with prefactor 4πr².
    """
    rhoatom::Vector{Float64}
    """
    Pseudo-atomic kinetic energy density on the radial grid.
    Note that there is no 4π nor r² prefactor here, **unlike rhoatom**.
    """
    tauatom::Union{Nothing,Vector{Float64}}
    "Spin-orbit coupling data, (ignored if `has_so` is false)"
    spin_orb::Union{Nothing,UpfSpinOrb}
    "PAW data, (ignored if `is_paw` is false)"
    paw::Union{Nothing,UpfPaw}
    "GIPAW data"
    gipaw::Union{Nothing,UpfGipaw}
    """
    Radii (bohr) where a table has a kink, by the name of its element (`PP_LOCAL`,
    `PP_NLCC`, `PP_TAUMOD`, `PP_RHOATOM`, `PP_TAUATOM`, `PP_BETA.i`, `PP_CHI.i`), from the
    attribute `breakpoints` of the element.
    An extension of UPF v2 (other readers ignore an unknown attribute): a generator that
    knows where its tables are only continuous (the matching radii, the radius of the local
    potential or of a model core) saves it, so that a reader can interpolate each side on
    its own. The writer also gives `breakpoint_indices`, the number of mesh points below
    each break. Empty if the file has none.
    """
    breakpoints::Dict{String,Vector{Float64}}
    UpfFile(fields...) = validated(new(fields...))
end

#= The specification: the UPF v2 documentation of QE, and what its reader (`upflib`)
requires. QE tolerates a header `mesh_size` that differs from the mesh (it takes the mesh),
and a negative occupation of a pseudo-wavefunction (not used for the starting
wavefunctions): those are only warnings (`@soft_check`). =#

const UPF_PSEUDO_TYPES = ("NC", "SL", "1/r", "US", "USPP", "PAW")
const UPF_RELATIVISTIC = ("no", "nonrelativistic", "scalar", "full")

function validate(h::UpfHeader)
    @check "UPF" haskey(PeriodicTable.elements, Symbol(h.element)) "element = $(repr(h.element)) is not a chemical symbol"
    @check "UPF" h.pseudo_type in UPF_PSEUDO_TYPES "pseudo_type = $(repr(h.pseudo_type)) is not one of $UPF_PSEUDO_TYPES"
    @check "UPF" isnothing(h.relativistic) || h.relativistic in UPF_RELATIVISTIC "relativistic = $(repr(h.relativistic)) is not one of $UPF_RELATIVISTIC"
    @check "UPF" h.is_paw == (h.pseudo_type == "PAW") "is_paw = $(h.is_paw), but pseudo_type = $(h.pseudo_type)"
    @check "UPF" h.is_ultrasoft == (h.pseudo_type in ("US", "USPP", "PAW")) "is_ultrasoft = $(h.is_ultrasoft), but pseudo_type = $(h.pseudo_type)"
    @check "UPF" !h.has_so || h.relativistic == "full" "has_so, but relativistic = $(repr(h.relativistic))"
    @check "UPF" h.z_valence > 0 "z_valence = $(h.z_valence)"
    @check "UPF" h.l_max >= -1 "l_max = $(h.l_max)"
    @check "UPF" h.mesh_size >= 1 && h.number_of_wfc >= 0 && h.number_of_proj >= 0 "mesh_size = $(h.mesh_size), number_of_wfc = $(h.number_of_wfc), number_of_proj = $(h.number_of_proj)"
    @check "UPF" something(h.wfc_cutoff, 0) >= 0 && something(h.rho_cutoff, 0) >= 0 "negative suggested cutoffs"
end

function validate(mesh::UpfMesh)
    n = length(mesh.r)
    @check "UPF" n >= 1 && mesh.r[1] >= 0 && all(i -> mesh.r[i] < mesh.r[i+1], 1:n-1) "PP_R is not increasing from r ≥ 0"
    @check "UPF" length(mesh.rab) == n "PP_RAB has $(length(mesh.rab)) values on a mesh of $n"
    @check "UPF" isnothing(mesh.mesh) || mesh.mesh == n "PP_R has $n values, but mesh = $(mesh.mesh)"
end

# i ≤ j, and the composite index ij = j(j-1)/2 + i (QE checks it)
function upf_check_indices(x, name)
    i, j, ij = x.first_index, x.second_index, x.composite_index
    isnothing(i) || isnothing(j) || @check "UPF" 1 <= i <= j "$name.$i.$j: indices not 1 ≤ i ≤ j"
    isnothing(i) || isnothing(j) || isnothing(ij) ||
        @check "UPF" ij == j * (j - 1) ÷ 2 + i "$name.$i.$j has composite_index $ij"
end
validate(x::UpfQij) = upf_check_indices(x, "PP_QIJ")
function validate(x::UpfQijl)
    @check "UPF" x.angular_momentum >= 0 "PP_QIJL: angular_momentum = $(x.angular_momentum)"
    upf_check_indices(x, "PP_QIJL")
end

function validate(aug::UpfAugmentation)
    @check "UPF" size(aug.q, 1) == size(aug.q, 2) "PP_Q is not square"
    @check "UPF" aug.q_with_l ? !isnothing(aug.qijls) : !isnothing(aug.qijs) "PP_AUGMENTATION: q_with_l = $(aug.q_with_l), but no PP_$(aug.q_with_l ? "QIJL" : "QIJ")"
    @check "UPF" aug.nqf >= 0 "nqf = $(aug.nqf)"
    aug.nqf > 0 &&
        @check "UPF" !isnothing(aug.qfcoefs) && !isnothing(aug.rinner) "nqf = $(aug.nqf), but no PP_QFCOEF or PP_RINNER"
    isnothing(aug.rinner) || isnothing(aug.nqlc) ||
        @check "UPF" length(aug.rinner) == aug.nqlc "PP_RINNER has $(length(aug.rinner)) values, but nqlc = $(aug.nqlc)"
end

validate(b::UpfBeta) = @check "UPF" b.angular_momentum >= 0 "PP_BETA: angular_momentum = $(b.angular_momentum)"

function validate(nl::UpfNonlocal)
    nb = length(nl.betas)
    @check "UPF" size(nl.dij) == (nb, nb) "PP_DIJ is $(size(nl.dij)) for $nb projectors"
    aug = nl.augmentation
    isnothing(aug) && return
    @check "UPF" size(aug.q) == (nb, nb) "PP_Q is $(size(aug.q)) for $nb projectors"
    for x in [something(aug.qijs, []); something(aug.qijls, [])]
        @check "UPF" something(x.second_index, 1) <= nb "an augmentation function of projector $(x.second_index), but $nb projectors"
    end
end

function validate(c::UpfChi)
    @check "UPF" c.l >= 0 "PP_CHI: l = $(c.l)"
    @soft_check "UPF" c.occupation >= 0 "PP_CHI: occupation = $(c.occupation) (QE: not used for the starting wavefunctions)"
    @check "UPF" isnothing(c.n) || c.n > c.l "PP_CHI: n = $(c.n), l = $(c.l)"
end

validate(x::UpfRelWfc) =
    @check "UPF" isnothing(x.lchi) || valid_j(x.lchi, x.jchi) "PP_RELWFC: jchi = $(x.jchi), lchi = $(x.lchi)"
validate(x::UpfRelBeta) =
    @check "UPF" isnothing(x.lll) || valid_j(x.lll, x.jjj) "PP_RELBETA: jjj = $(x.jjj), lll = $(x.lll)"

validate(w::UpfWfc) = @check "UPF" w.l >= 0 && w.index >= 1 "a full wavefunction with l = $(w.l), index = $(w.index)"

validate(o::UpfGipawCoreOrbital) =
    @check "UPF" 0 <= o.l < o.n "PP_GIPAW_CORE_ORBITAL: n = $(o.n), l = $(o.l)"

# The names of the tables of a file (see `breakpoints`)
function upf_table_names(file::UpfFile)
    names = [name for (name, x) in ("PP_LOCAL" => file.local_, "PP_NLCC" => file.nlcc,
                                    "PP_TAUMOD" => file.taumod, "PP_RHOATOM" => file.rhoatom,
                                    "PP_TAUATOM" => file.tauatom, "PP_DIJ" => file.nonlocal.dij)
             if !isnothing(x)]
    indexed(name, i) = isnothing(i) ? name : "$name.$i"
    append!(names, indexed("PP_BETA", b.index) for b in file.nonlocal.betas)
    append!(names, indexed("PP_CHI", c.index) for c in file.pswfc)
    return names
end

function validate(file::UpfFile)
    h = file.header
    n = length(file.mesh.r)
    nb = length(file.nonlocal.betas)
    on_mesh(x) = isnothing(x) || length(x) == n
    for (name, table) in ("PP_LOCAL" => file.local_, "PP_NLCC" => file.nlcc,
                          "PP_TAUMOD" => file.taumod, "PP_RHOATOM" => file.rhoatom,
                          "PP_TAUATOM" => file.tauatom)
        @check "UPF" on_mesh(table) "$name has $(length(table)) values on a mesh of $n"
    end
    @check "UPF" all(b -> on_mesh(b.beta), file.nonlocal.betas) "a PP_BETA is not on the mesh"
    @check "UPF" all(c -> on_mesh(c.chi), file.pswfc) "a PP_CHI is not on the mesh"
    @check "UPF" h.is_coulomb == isnothing(file.local_) "is_coulomb = $(h.is_coulomb), but PP_LOCAL is $(isnothing(file.local_) ? "missing" : "present")"
    @check "UPF" h.core_correction == !isnothing(file.nlcc) "core_correction = $(h.core_correction), but PP_NLCC is $(isnothing(file.nlcc) ? "missing" : "present")"
    @check "UPF" h.with_metagga_info == !isnothing(file.taumod) == !isnothing(file.tauatom) "with_metagga_info = $(h.with_metagga_info), but PP_TAUMOD and PP_TAUATOM are not both $(h.with_metagga_info ? "present" : "missing")"
    @check "UPF" h.has_so == !isnothing(file.spin_orb) "has_so = $(h.has_so), but PP_SPIN_ORB is $(isnothing(file.spin_orb) ? "missing" : "present")"
    @check "UPF" h.is_paw == !isnothing(file.paw) "is_paw = $(h.is_paw), but PP_PAW is $(isnothing(file.paw) ? "missing" : "present")"
    @check "UPF" h.has_gipaw == !isnothing(file.gipaw) "has_gipaw = $(h.has_gipaw), but PP_GIPAW is $(isnothing(file.gipaw) ? "missing" : "present")"
    @check "UPF" h.has_wfc == !isnothing(file.full_wfc) "has_wfc = $(h.has_wfc), but PP_FULL_WFC is $(isnothing(file.full_wfc) ? "missing" : "present")"
    @check "UPF" (h.is_ultrasoft || h.is_paw) == !isnothing(file.nonlocal.augmentation) "pseudo_type = $(h.pseudo_type), but PP_AUGMENTATION is $(isnothing(file.nonlocal.augmentation) ? "missing" : "present")"
    @check "UPF" h.number_of_proj == nb "number_of_proj = $(h.number_of_proj), but $nb PP_BETA"
    @soft_check "UPF" h.mesh_size == n "mesh_size = $(h.mesh_size), but a mesh of $n (QE takes the mesh)"
    @check "UPF" h.number_of_wfc == length(file.pswfc) "number_of_wfc = $(h.number_of_wfc), but $(length(file.pswfc)) PP_CHI"
    @check "UPF" all(b -> b.angular_momentum <= h.l_max, file.nonlocal.betas) "a projector of l > l_max = $(h.l_max)"
    aug = file.nonlocal.augmentation
    if !isnothing(aug)
        @check "UPF" all(x -> on_mesh(x.qij), something(aug.qijs, [])) "a PP_QIJ is not on the mesh"
        @check "UPF" all(x -> on_mesh(x.qijl), something(aug.qijls, [])) "a PP_QIJL is not on the mesh"
    end
    so = file.spin_orb
    if !isnothing(so)
        @check "UPF" length(so.relwfcs) == length(file.pswfc) "$(length(so.relwfcs)) PP_RELWFC for $(length(file.pswfc)) PP_CHI"
        @check "UPF" length(so.relbetas) == nb "$(length(so.relbetas)) PP_RELBETA for $nb PP_BETA"
        @check "UPF" all(((r, b),) -> isnothing(r.lll) || r.lll == b.angular_momentum, zip(so.relbetas, file.nonlocal.betas)) "the l of a PP_RELBETA is not that of its PP_BETA"
    end
    paw = file.paw
    if !isnothing(paw)
        @check "UPF" length(paw.occupations) == nb "PP_OCCUPATIONS has $(length(paw.occupations)) values for $nb projectors"
        @check "UPF" on_mesh(paw.ae_nlcc) && on_mesh(paw.ae_vloc) "PP_AE_NLCC or PP_AE_VLOC is not on the mesh"
        @check "UPF" length(paw.aewfcs) in (0, nb) && length(paw.pswfcs) in (0, nb) "PP_PAW: not one wavefunction per projector"
    end
    fw = file.full_wfc
    if !isnothing(fw)
        @check "UPF" length(fw.aewfcs) == nb && length(fw.pswfcs) in (0, nb) "PP_FULL_WFC: not one wavefunction per projector"
        @check "UPF" all(w -> on_mesh(w.wfc), [fw.aewfcs; fw.pswfcs]) "a PP_AEWFC or PP_PSWFC is not on the mesh"
    end
    isnothing(file.gipaw) ||
        @check "UPF" all(o -> on_mesh(o.core_orbital), file.gipaw.core_orbitals) "a PP_GIPAW_CORE_ORBITAL is not on the mesh"
    tables = upf_table_names(file)
    for (name, breaks) in file.breakpoints
        @check "UPF" name in tables "breakpoints of $name, which the file does not have"
        @check "UPF" all(b -> file.mesh.r[1] < b < file.mesh.r[end], breaks) "a breakpoint of $name is not inside the mesh"
        @soft_check "UPF" issorted(breaks) "the breakpoints of $name are not increasing"
    end
end

# the fields of UPF proper; no breakpoints
UpfFile(identifier, version, info, header, mesh, nlcc, taumod, local_, nonlocal, pswfc,
        full_wfc, rhoatom, tauatom, spin_orb, paw, gipaw) =
    UpfFile(identifier, version, info, header, mesh, nlcc, taumod, local_, nonlocal, pswfc,
            full_wfc, rhoatom, tauatom, spin_orb, paw, gipaw, Dict{String,Vector{Float64}}())

function UpfFile(path::AbstractString; identifier="", validation=VALIDATION[])
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return UpfFile(io; identifier, validation)
    end
end

function UpfFile(io::IO; identifier="", validation=VALIDATION[])
    version = _get_upf_version(io)
    with_validation(validation; source=identifier) do
        version == 1 && return _from_dict(UpfFile, upf1_to_dict(read(io, String)); identifier)
        return upf2_parse_psp(io; identifier)
    end
end

"""
Write UPF file to `io`. Note that the file is written the UPF version passed
as a string and this value takes preference over the version stored in `psp`.

Thread-safe.
"""
function save_psp(io::IO, psp::UpfFile; format="UPF v2.0.1")
    if format == "UPF v2.0.1"
        upf2_dump_psp(io, psp)
    elseif format == "UPF v1.old"
        error("$format not supported yet.")
    else
        error("$format is unknown. Only 'UPF v2.0.1' and 'UPF v1.old' are known.")
    end
end

function _get_upf_version(io::IO)::Int
    pos = position(io)
    line1 = readline(io)
    line2 = readline(io)
    seek(io, pos)
    if occursin("<PP_INFO>", line1)
        # Old UPF files start with the `<PP_INFO>` section
        return 1
    elseif occursin("UPF version=\"2.0.1\"", line1)
        # New UPF files with schema are in XML and start with a version tag
        return 2
    elseif occursin("xml version=\"1.0\"", line1) && occursin("UPF version=\"2.0.1\"", line2)
        # dumped UPF files start with a xml declaration
        return 2
    else
        error("Unknown UPF version")
    end
end

function _get_upf_version(path::AbstractString)::Int
    open(path, "r") do io
        return _get_upf_version(io)
    end
end

identifier(psp::UpfFile)::String = psp.identifier
format(file::UpfFile)::String = "UPF v$(file.version)"
"""
Libxc string of a QE `XC-000i-000i-...` functional (inverse of `libxc_to_qe_libxc`). The slots
are LDA exchange, LDA correlation, GGA exchange, GGA correlation, meta-GGA exchange, meta-GGA
correlation; `L` marks a Libxc id, `I` a QE-internal index (zero: no term). `"UNKNOWN"` for a
QE-internal meta-GGA or an id without a Libxc name.
"""
function qe_libxc_to_libxc(functional::AbstractString)::String
    slot_tables = (UPF_LDA_EXCHANGE, UPF_LDA_CORRELATION, UPF_GGA_EXCHANGE,
                   UPF_GGA_CORRELATION)
    terms = String[]
    for (slot, term) in enumerate(split(uppercase(strip(functional)), '-')[2:end])
        m = match(r"^(\d+)([LI])$", term)
        (isnothing(m) || slot > 6) && return "UNKNOWN"
        id = parse(Int, m[1])
        iszero(id) && continue
        if m[2] == "L"
            name = get(LIBXC_FUNCTIONALS_BY_ID, id, nothing)
        else
            slot > 4 && return "UNKNOWN"
            entry = findfirst(e -> e["i"] == id, slot_tables[slot])
            name = isnothing(entry) ? nothing : slot_tables[slot][entry]["libxc"]
        end
        (isnothing(name) || isempty(name)) && return "UNKNOWN"
        push!(terms, name)
    end
    return isempty(terms) ? "UNKNOWN" : join(terms, ' ')
end

function libxc_string(header::UpfHeader)
    functional = header.functional
    startswith(uppercase(strip(functional)), "XC-") && return qe_libxc_to_libxc(functional)
    upf_codes = lowercase.(split(functional))

    if length(upf_codes) == 1  # Short code
        # special case r2SCAN01 as produced by METAPSP,
        # since QE so far only expects LDA and GGA XCs
        upf_codes[1] == "r2scan01" && return "mgga_x_r2scan01 mgga_c_r2scan01"

        entry_index = findfirst(e -> e["name"] == upf_codes[1], UPF_SHORT_NAMES)
        isnothing(entry_index) && return "UNKNOWN"
        long_code = UPF_SHORT_NAMES[entry_index]["full_name"]
        upf_codes = split(long_code, '+')
    end

    @assert length(upf_codes) == 4

    if upf_codes[3] == "nogx"
        exc_index = findfirst(e -> e["name"] == upf_codes[1], UPF_LDA_EXCHANGE)
        exc = UPF_LDA_EXCHANGE[exc_index]["libxc"]
    else
        exc_index = findfirst(e -> e["name"] == upf_codes[3], UPF_GGA_EXCHANGE)
        exc = UPF_GGA_EXCHANGE[exc_index]["libxc"]
    end

    if upf_codes[4] == "nogc"
        corr_index = findfirst(e -> e["name"] == upf_codes[2], UPF_LDA_CORRELATION)
        corr = UPF_LDA_CORRELATION[corr_index]["libxc"]
    else
        corr_index = findfirst(e -> e["name"] == upf_codes[4], UPF_GGA_CORRELATION)
        corr = UPF_GGA_CORRELATION[corr_index]["libxc"]
    end

    return join([exc, corr], ' ')
end
functional(file::UpfFile)::Vector{Symbol} = map(s -> Symbol(lowercase(s)), split(libxc_string(file.header)))
element(file::UpfFile)::PeriodicTable.Element = PeriodicTable.elements[Symbol(file.header.element)]
is_norm_conserving(file::UpfFile)::Bool = file.header.pseudo_type == "NC"
is_ultrasoft(file::UpfFile)::Bool = file.header.pseudo_type in ("US", "USPP")
is_paw(file::UpfFile)::Bool = file.header.pseudo_type == "PAW"
has_spin_orbit(file::UpfFile)::Bool = file.header.has_so
has_model_core_charge_density(file::UpfFile)::Bool = file.header.core_correction
valence_charge(file::UpfFile)::Float64 = file.header.z_valence

function libxc_to_qe_libxc(libxc_string::AbstractString)::String
    # https://gitlab.com/QEF/q-e/-/blob/develop/Modules/funct.f90?ref_type=heads#L301
    # ! NOTE FOR LIBXC USERS: to use libxc functionals you must enforce them from input (use
    # ! 'input_dft' in &system) and write their IDs in the input string. The only notation
    # ! now allowed (v7.0) for input DFTs containing Libxc terms is:
    # ! XC-000i-000i-000i-000i-000i-000i
    # ! where you put the functional IDs instead of the zeros and an 'L' instead of
    # ! 'i' if the functional is from Libxc. The order is the usual one:
    # ! LDAexch - LDAcorr - GGAexch - GGAcorr - MGGAexch - MGGAcorr
    # ! however QE will automatically adjust it if needed. You can skip zero tails (e.g.
    # ! you don't need GGA/MGGA slots if the dft is LDA only and so on.
    # ! You can use combinations of qe and libxc functionals, when they are compatible.
    # ! You can also add vdW terms after it, for example, sla+pw+rw86+vdw2 is:
    # ! input_dft='XC-001i-004i-013i-vdw2'.
    # ! For more details see the user_guide (in 'Doc' folder).
    # https://www.quantum-espresso.org/wp-content/uploads/2022/03/user_guide.pdf
    # Some functionals in libxc incorporate the exchange part and the correlation one into one term
    # only (e.g. the ones that include the ‘ xc’ kind label in their name). In these cases the whole
    # functional is formally treated as ‘correlation only’ by Quantum ESPRESSO. This does not
    # imply any loss of information.
    terms = zeros(Int, 6)
    for substring in split(libxc_string)
        substring_parts = split(substring, '_')
        if substring_parts[1] == "LDA"
            if substring_parts[2] == "X"
                terms[1] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[2] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        elseif substring_parts[1] == "GGA"
            if substring_parts[2] == "X"
                terms[3] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[4] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        elseif substring_parts[1] == "MGGA"
            if substring_parts[2] == "X"
                terms[5] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            elseif substring_parts[2] in ("C", "XC")
                terms[6] = LIBXC_FUNCTIONALS_BY_NAME[substring]
            else
                error("Unsupported LibXC functional $(libxc_string)")
            end
        else
            error("Unsupported LibXC functional $(libxc_string)")
        end
    end
    # zero tails can be skipped
    n_terms = something(findlast(!iszero, terms), 0)
    return "XC" * prod(terms[1:n_terms]; init="") do term
        iszero(term) && return @sprintf "-%03dI" term
        return @sprintf "-%03dL" term
    end
end

function libxc_to_qe_internal(libxc_string::AbstractString, sep::Char=' ')::String
    luts = Dict(
        "LDA" => Dict("X" => UPF_LDA_EXCHANGE, "C" => UPF_LDA_CORRELATION),
        "GGA" => Dict("X" => UPF_GGA_EXCHANGE, "C" => UPF_GGA_CORRELATION),
    )
    labels_indices = Dict(
        "LDA" => Dict("X" => 1, "C" => 2),
        "GGA" => Dict("X" => 3, "C" => 4)
    )
    qe_labels = ["nox", "noc", "nogx", "nogc"]
    for libxc_name in split(libxc_string)
        family, x_or_c = split(libxc_name, '_')[1:2]

        family_lut = get(luts, family, nothing)
        isnothing(family_lut) && error("Could not convert $(libxc_string) to QE complete name.")

        lut = get(family_lut, x_or_c, nothing)
        isnothing(lut) && return error("Could not convert $(libxc_string) to QE complete name.")

        labels_idx = labels_indices[family][x_or_c]

        lut_idx = findfirst(entry -> entry["libxc"] == libxc_name, lut)
        qe_labels[labels_idx] = isnothing(lut_idx) ? nothing : lut[lut_idx]["name"]
    end
    # LibXC GGA functionals include LDA exchange and correlation: Slater exchange, Perdew-Wang correlation.
    # Lee-Yang-Parr functionals are an exception, but they don't have QE old equivalents in the lookup tables,
    # so we don't handle them here.
    if qe_labels[3] != "nogx"

        qe_labels[1] = "sla"
    end
    if qe_labels[4] != "nogc"
        qe_labels[2] = "pw"
    end
    any(isnothing.(qe_labels)) && error("Could not convert $(libxc_string) to QE complete name.")
    return join(qe_labels, sep)
end
