const MAX_N_BETAS = 20
const MAX_N_CHIS = 20
const MAX_N_FULL_WFCS = 20
const MAX_N_GIPAW_CORE_ORBITALS = 20
const MAX_N_GIPAW_ORBITALS = 20

abstract type AbstractPpVector end

struct PpVector <: AbstractPpVector
    _::Vector{Float64}
end

function Serde.deser(::Type{<:AbstractPpVector}, ::Type{Vector{Float64}}, text::String)
    return parse.(Float64, strip.(split(text)))
end

struct PpInputfile
    _::String
end

# Note: PP_INFO usually has some text before the PP_INPUTFILE tag, which is lost
struct PpInfo
    _::Union{Nothing,String}
    PP_INPUTFILE::Union{Nothing,PpInputfile}
end

struct PpHeader
    generated::Union{Nothing,String}
    author::Union{Nothing,String}
    date::Union{Nothing,String}
    comment::Union{Nothing,String}
    element::String
    pseudo_type::String
    relativistic::Union{Nothing,String}
    is_ultrasoft::Bool
    is_paw::Bool
    is_coulomb::Bool
    has_so::Bool
    has_wfc::Bool
    has_gipaw::Bool
    paw_as_gipaw::Union{Nothing,Bool}
    core_correction::Bool
    functional::String
    z_valence::Float64
    total_psenergy::Union{Nothing,Float64}
    wfc_cutoff::Union{Nothing,Float64}
    rho_cutoff::Union{Nothing,Float64}
    l_max::Int
    l_max_rho::Union{Nothing,Int}
    l_local::Union{Nothing,Int}
    mesh_size::Int
    number_of_wfc::Int
    number_of_proj::Int
end

Serde.deser(::Type{PpHeader}, ::Type{Bool}, text::String) = occursin("t", lowercase(text))

struct PpMesh
    dx::Union{Nothing,Float64}
    xmin::Union{Nothing,Float64}
    rmax::Union{Nothing,Float64}
    zmesh::Union{Nothing,Float64}
    PP_R::PpVector
    PP_RAB::PpVector
end

struct PpBeta <: AbstractPpVector
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    angular_momentum::Union{Nothing,Int}
    cutoff_radius_index::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    _::Vector{Float64}
end

struct PpDij
    rows::Union{Nothing,Int}
    columns::Union{Nothing,Int}
    _::Matrix{Float64}
end

function Serde.deser(::Type{PpDij}, ::Type{Matrix{Float64}}, text::String)
    vector = parse.(Float64, strip.(split(text)))
    if sqrt(length(vector)) != floor(sqrt(length(vector)))
        throw(ArgumentError("PP_DIJ should be a square matrix"))
    end
    n = Int(sqrt(length(vector)))
    return reshape(vector, n, n)
end

struct PpQij <: AbstractPpVector
    first_index::Int
    second_index::Int
    composite_index::Int
    _::Vector{Float64}
end

@eval struct PpAugmentation
    q_with_l::Bool
    nqf::Int
    nqlc::Int
    PP_Q::PpVector
    $([:($(Symbol("PP_QIJ_", i, "_", j))::Union{Nothing,PpQij})
       for i in 1:MAX_N_BETAS for j in i:MAX_N_BETAS]...)
end
for i in 1:MAX_N_BETAS
    for j in 1:MAX_N_BETAS
        @eval function Serde.custom_name(
                ::Type{PpAugmentation},
                ::Val{$Symbol("PP_QIJ_", $i, "_", $j)}
            )
            string("PP_QIJ.", $i, ".", $j)
        end
    end
end

Serde.deser(::Type{PpAugmentation}, ::Type{Bool}, t::String) = occursin("t", lowercase(t))

# Builds a struct with fields PP_BETA1, PP_BETA2, ..., PP_BETAN, where N = MAX_N_BETAS
# and PP_DIJ
@eval struct PpNonlocal
    $([:($(Symbol("PP_BETA_", i))::Union{Nothing,PpBeta})
       for i in 1:MAX_N_BETAS]...)
    PP_DIJ::PpDij
    PP_AUGMENTATION::Union{Nothing,PpAugmentation}
end
# Tells Serde that in the XML file the field PP_BETA$(i) is named "PP_BETA.$(i)"
for i in 1:MAX_N_BETAS
    @eval function Serde.custom_name(
            ::Type{PpNonlocal},
            ::Val{$Symbol("PP_BETA_", $i)}
        )
        string("PP_BETA.", $i)
    end
end

struct PpChi <: AbstractPpVector
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    l::Union{Nothing,Int}
    occupation::Union{Nothing,Float64}
    n::Union{Nothing,Int}
    cutoff_radius::Union{Nothing,Float64}
    ultrasoft_cutoff_radius::Union{Nothing,Float64}
    _::Vector{Float64}
end

# Same as above for PP_NONLOCAL and PP_BETA but for PP_PSWFC and PP_CHI
@eval struct PpPswfc
    $([:($(Symbol("PP_CHI_", i))::Union{Nothing,PpChi})
       for i in 1:MAX_N_CHIS]...)
end
for i in 1:MAX_N_CHIS
    @eval function Serde.custom_name(
            ::Type{PpPswfc},
            ::Val{$Symbol("PP_CHI_", $i)}
        )
        string("PP_CHI.", $i)
    end
end

struct PpFullWfcAewfc <: AbstractPpVector
    l::Union{Nothing,Int}
    label::Union{Nothing,String}
    _::Vector{Float64}
end

struct PpFullWfcPswfc <: AbstractPpVector
    l::Union{Nothing,Int}
    label::Union{Nothing,String}
    _::Vector{Float64}
end

@eval struct PpFullWfc
    number_of_wfc::Int
    $([
        [:($(Symbol("PP_AEWFC_", i))::Union{Nothing,PpFullWfcAewfc})
         for i in 1:MAX_N_FULL_WFCS];
        [:($(Symbol("PP_PSWFC_", i))::Union{Nothing,PpFullWfcPswfc})
         for i in 1:MAX_N_FULL_WFCS]
    ]...)
end
for i in 1:MAX_N_FULL_WFCS
    @eval function Serde.custom_name(
            ::Type{PpFullWfc},
            ::Val{$Symbol("PP_AEWFC_", $i)}
        )
        string("PP_AEWFC.", $i)
    end
    @eval function Serde.custom_name(
            ::Type{PpFullWfc},
            ::Val{$Symbol("PP_PSWFC_", $i)}
        )
        string("PP_PSWFC.", $i)
    end
end

struct PpPaw
    paw_data_format::Int
    core_energy::Union{Nothing,Float64}
    PP_OCCUPATIONS::PpVector
    PP_AE_NLCC::PpVector
    PP_AE_VLOC::PpVector
end

struct PpRelbeta
    index::Int
    lll::Int
    jjj::Float64
end

struct PpRelwfc
    index::Int
    lchi::Int
    jchi::Float64
    nn::Int
end

@eval struct PpSpinOrb
    $([
        [:($(Symbol("PP_RELBETA_", i))::Union{Nothing,PpRelbeta})
         for i in 1:MAX_N_BETAS];
        [:($(Symbol("PP_RELWFC_", i))::Union{Nothing,PpRelwfc})
         for i in 1:MAX_N_CHIS]
    ]...)
end
for i in 1:MAX_N_BETAS
    @eval function Serde.custom_name(
            ::Type{PpSpinOrb},
            ::Val{$Symbol("PP_RELBETA_", $i)}
        )
        string("PP_RELBETA.", $i)
    end
end
for i in 1:MAX_N_CHIS
    @eval function Serde.custom_name(
            ::Type{PpSpinOrb},
            ::Val{$Symbol("PP_RELWFC_", $i)}
        )
        string("PP_RELWFC.", $i)
    end
end

struct PpGipawCoreOrbital <: AbstractPpVector
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    n::Union{Nothing,Int}
    l::Union{Nothing,Int}
    _::Vector{Float64}
end

Serde.deser(::Type{PpGipawCoreOrbital}, ::Type{Int}, t::String) = Int(parse(Float64, t))

# Same as above for PP_NONLOCAL, PP_PSWFC
@eval struct PpGipawCoreOrbitals
    $([:($(Symbol("PP_GIPAW_CORE_ORBITAL_", i))::Union{Nothing,PpGipawCoreOrbital})
       for i in 1:MAX_N_GIPAW_CORE_ORBITALS]...)
end
for i in 1:MAX_N_GIPAW_CORE_ORBITALS
    @eval function Serde.custom_name(
            ::Type{PpGipawCoreOrbitals},
            ::Val{$Symbol("PP_GIPAW_CORE_ORBITAL_", $i)}
        )
        string("PP_GIPAW_CORE_ORBITAL.", $i)
    end
end

struct PpGipawOrbital <: AbstractPpVector
    index::Union{Nothing,Int}
    label::Union{Nothing,String}
    l::Union{Nothing,Int}
    PP_GIPAW_WFS_AE::PpVector
    PP_GIPAW_WFS_PS::PpVector
end

# Same as above for PP_NONLOCAL, PP_PSWFC, and PP_GIPAW_CORE_ORBITALS
@eval struct PpGipawOrbitals
    $([:($(Symbol("PP_GIPAW_ORBITAL_", i))::Union{Nothing,PpGipawOrbital})
       for i in 1:MAX_N_GIPAW_ORBITALS]...)
end
for i in 1:MAX_N_GIPAW_ORBITALS
    @eval function Serde.custom_name(
            ::Type{PpGipawOrbitals},
            ::Val{$Symbol("PP_GIPAW_ORBITAL_", $i)}
        )
        string("PP_GIPAW_ORBITAL.", $i)
    end
end

struct PpGipawVlocal
    PP_GIPAW_VLOCAL_AE::PpVector
    PP_GIPAW_VLOCAL_PS::PpVector
end

struct PpGipaw
    gipaw_data_format::Int
    PP_GIPAW_CORE_ORBITALS::PpGipawCoreOrbitals
    PP_GIPAW_ORBITALS::Union{Nothing,PpGipawOrbitals}
    PP_GIPAW_VLOCAL::Union{Nothing,PpGipawVlocal}
end

struct UPFSerde
    version::String
    PP_INFO::PpInfo
    PP_HEADER::PpHeader
    PP_MESH::PpMesh
    PP_NLCC::Union{Nothing,PpVector}
    PP_LOCAL::Union{Nothing,PpVector}
    PP_NONLOCAL::PpNonlocal
    PP_PSWFC::Union{Nothing,PpPswfc}
    PP_FULL_WFC::Union{Nothing,PpFullWfc}
    PP_RHOATOM::Union{Nothing,PpVector}
    PP_SPIN_ORB::Union{Nothing,PpSpinOrb}
    PP_PAW::Union{Nothing,PpPaw}
    PP_GIPAW::Union{Nothing,PpGipaw}
end
