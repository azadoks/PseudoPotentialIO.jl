"""
ABINIT PSeudoPotential format 8 header block.
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
    "Signals presence of format extensions: 0 (none), 1 (valence charge), 2 (spin-orbit), 3 (spin-orbit and valence charge)"
    extension_switch::Int
    "Number of spin-orbit projectors for each angular momentum (if present; 0 for l = 0)"
    nprojso::Union{Nothing,Vector{Int}}
    Psp8Header(fields...) = validated(new(fields...))
end

function validate(h::Psp8Header)
    @check "PSP8" h.pspcod == 8 "pspcod = $(h.pspcod)"
    @check "PSP8" 0 < h.zion <= h.zatom "zion = $(h.zion), zatom = $(h.zatom)"
    @check "PSP8" psp8_known_pspxc(h.pspxc) "pspxc = $(h.pspxc) is not a known functional"
    @check "PSP8" h.lmax >= 0 && h.lloc >= 0 "lmax = $(h.lmax), lloc = $(h.lloc)"
    @check "PSP8" h.mmax >= 1 "mmax = $(h.mmax)"
    @check "PSP8" h.rchrg >= 0 && h.fchrg >= 0 "rchrg = $(h.rchrg), fchrg = $(h.fchrg)"
    @check "PSP8" length(h.nproj) == h.lmax + 1 && all(>=(0), h.nproj) "nproj = $(h.nproj) for lmax = $(h.lmax)"
    @check "PSP8" h.extension_switch in 0:3 "extension_switch = $(h.extension_switch)"
    has_so = h.extension_switch in (2, 3)
    @check "PSP8" isnothing(h.nprojso) == !has_so "extension_switch = $(h.extension_switch), but nprojso = $(h.nprojso)"
    isnothing(h.nprojso) ||
        @check "PSP8" length(h.nprojso) == h.lmax + 1 && h.nprojso[1] == 0 && all(>=(0), h.nprojso) "nprojso = $(h.nprojso) for lmax = $(h.lmax)"
end

function Base.show(io::IO, ::MIME"text/plain", header::Psp8Header)
    println(io, typeof(header))
    @printf "%032s: %s\n" "description" header.title
    @printf "%032s: %0.6f\n" "atomic charge" header.zatom
    @printf "%032s: %0.6f\n" "valence charge" header.zion
    @printf "%032s: %06d\n" "generation date" header.pspd
    @printf "%032s: %d\n" "PSP format version" header.pspcod
    if header.pspxc > 0
        @printf "%032s: %d\n" "exchange-correlation (ABINIT)" header.pspxc
    else
        @printf "%032s: %d\n" "exchange-correlation (LibXC)" header.pspxc
    end
    @printf "%032s: %d\n" "max angular momentum" header.lmax
    @printf "%032s: %d\n" "local potential angular momentum" header.lloc
    @printf "%032s: %d\n" "radial mesh size" header.mmax
    @printf "%032s: %0.6f\n" "r2well" header.r2well
    @printf "%032s: %0.6f\n" "core charge radius" header.rchrg
    @printf "%032s: %0.6f\n" "model core charge present" header.fchrg
    @printf "%032s: %0.6f\n" "qchrg" header.qchrg
    @printf "%032s: %s\n" "number of projectors" join(header.nproj, ", ")
    @printf "%032s: %d\n" "extension switch" header.extension_switch
    if header.extension_switch in (2, 3)
        @printf "%032s: %s\n" "number of spin-orbit projectors" join(header.nprojso, ", ")
    end
end

"""
ABINIT PSeudoPotential format 8 file contents. Information on the file format specification
and the meaning of the quantities within the file can be found on the
["psp8" page](https://docs.abinit.org/developers/psp8_info/) of the ABINIT documentation.
"""
struct Psp8File <: PsPFile
    "Identifier"
    identifier::String
    "Various pseudopotential metadata"
    header::Psp8Header
    "Uniform radial grid starting at `r = 0.0`"
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
    "Model core charge density with 4π prefactor"
    rhoc::Union{Nothing,Vector{Float64}}
    "First derivative of the model core charge density"
    d_rhoc_dr::Union{Nothing,Vector{Float64}}
    "Second derivative of the model core charge density"
    d2_rhoc_dr2::Union{Nothing,Vector{Float64}}
    "Third derivative of the model core charge density"
    d3_rhoc_dr3::Union{Nothing,Vector{Float64}}
    "Fourth derivative of the model core charge density"
    d4_rhoc_dr4::Union{Nothing,Vector{Float64}}
    "Valence charge density with 4π prefactor"
    rhov::Union{Nothing,Vector{Float64}}
    "All-electron valence charge density with 4π prefactor"
    ae_rhov::Union{Nothing,Vector{Float64}}
    "All-electron core charge density with 4π prefactor"
    ae_rhoc::Union{Nothing,Vector{Float64}}
    Psp8File(fields...) = validated(new(fields...))
end

# projectors and energies by l, as the header says
function psp8_check_projectors(projectors, ekb, nproj, mmax, name)
    @check "PSP8" length(projectors) == length(ekb) == length(nproj) "$name: $(length(projectors)) l, but $(length(nproj)) in the header"
    for (l, (p, e, n)) in enumerate(zip(projectors, ekb, nproj))
        @check "PSP8" length(p) == length(e) == n "$name: $(length(p)) projectors and $(length(e)) energies for l = $(l - 1), but $n in the header"
        @check "PSP8" all(x -> length(x) == mmax, p) "$name: a projector of l = $(l - 1) is not on the $mmax points of the mesh"
    end
end

function validate(file::Psp8File)
    h = file.header
    mmax = h.mmax
    on_mesh(x) = isnothing(x) || length(x) == mmax
    @check "PSP8" length(file.rgrid) == mmax "$(length(file.rgrid)) mesh points, but mmax = $mmax"
    # uniform, from 0
    step = length(file.rgrid) > 1 ? file.rgrid[2] - file.rgrid[1] : 0.0
    @check "PSP8" !isempty(file.rgrid) && iszero(file.rgrid[1]) && step > 0 &&
                  all(i -> abs(file.rgrid[i] - (i - 1) * step) <= 1e-8 * file.rgrid[end], eachindex(file.rgrid)) "the radial mesh is not uniform from r = 0"
    @check "PSP8" on_mesh(file.v_local) "the local potential is not on the mesh"
    psp8_check_projectors(file.projectors, file.ekb, h.nproj, mmax, "projectors")
    @check "PSP8" isnothing(file.projectors_so) == isnothing(file.ekb_so) == isnothing(h.nprojso) "spin-orbit projectors, but not their energies or their number (or the reverse)"
    isnothing(file.projectors_so) || isnothing(h.nprojso) ||
        psp8_check_projectors(file.projectors_so, file.ekb_so, h.nprojso, mmax, "spin-orbit projectors")
    nlcc = [getfield(file, Symbol(k)) for k in PSP8_NLCC_COLUMNS]
    @check "PSP8" all(isnothing, nlcc) == !(h.fchrg > 0) && all(on_mesh, nlcc) "fchrg = $(h.fchrg), but the model core charge has $(count(!isnothing, nlcc)) of its 5 columns on the mesh"
    @check "PSP8" isnothing(file.rhov) == !(h.extension_switch in (1, 3)) "extension_switch = $(h.extension_switch), but the valence density is $(isnothing(file.rhov) ? "missing" : "present")"
    @check "PSP8" !isnothing(file.rhov) || isnothing(file.ae_rhov) && isnothing(file.ae_rhoc) "all-electron densities without the valence density"
    @check "PSP8" on_mesh(file.rhov) && on_mesh(file.ae_rhov) && on_mesh(file.ae_rhoc) "a valence density is not on the mesh"
end

Psp8File(io::IO; identifier="", validation=VALIDATION[]) =
    with_validation(validation; source=identifier) do
        return _from_dict(Psp8File, psp8_to_dict(read(io, String)); identifier)
    end

function Psp8File(path::AbstractString; identifier="", validation=VALIDATION[])
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return Psp8File(io; identifier, validation)
    end
end

# The Libxc ids of a negative `pspxc`: its digits, three by three (`-001012` is 1 and 12)
function psp8_libxc_ids(pspxc::Integer)
    libxc_digits = digits(-pspxc)  # least significant first
    n_ids = cld(length(libxc_digits), 3)
    append!(libxc_digits, zeros(Int, 3n_ids - length(libxc_digits)))
    return [evalpoly(10, libxc_digits[3i-2:3i]) for i in n_ids:-1:1]
end

psp8_known_pspxc(pspxc::Integer) =
    pspxc > 0 ? any(entry -> entry["i"] == pspxc, PSP8_EXCHANGE_CORRELATION) :
    pspxc < 0 && all(in(keys(LIBXC_FUNCTIONALS_BY_ID)), psp8_libxc_ids(pspxc))

# The Libxc names of an ABINIT `pspxc` (also of the HGH formats of ABINIT); empty if unknown
libxc_string(header::Psp8Header)::String = abinit_libxc_string(header.pspxc)

function abinit_libxc_string(pspxc::Integer)::String
    psp8_known_pspxc(pspxc) || return ""
    if pspxc > 0  # ABINIT code
        entry_index = findfirst(entry -> entry["i"] == pspxc, PSP8_EXCHANGE_CORRELATION)
        return PSP8_EXCHANGE_CORRELATION[entry_index]["libxc"]
    end
    return join((LIBXC_FUNCTIONALS_BY_ID[id] for id in psp8_libxc_ids(pspxc)), ' ')
end

identifier(file::Psp8File)::String = file.identifier
format(::Psp8File)::String = "PSP8"
element(file::Psp8File)::PeriodicTable.Element = PeriodicTable.elements[Int(file.header.zatom)]
functional(file::Psp8File)::Vector{Symbol} = map(s -> Symbol(lowercase(s)), split(libxc_string(file.header)))
valence_charge(file::Psp8File)::Float64 = file.header.zion
is_norm_conserving(file::Psp8File)::Bool = true
is_ultrasoft(file::Psp8File)::Bool = false
is_paw(file::Psp8File)::Bool = false
has_spin_orbit(file::Psp8File)::Bool = file.header.extension_switch in (2, 3)
has_model_core_charge_density(file::Psp8File)::Bool = !isnothing(file.rhoc)

function libxc_to_abinit_libxc(libxc_string::AbstractString)::Int
    abinit_libxc_int_str = "-" * prod(split(libxc_string)) do substring
                                 return @sprintf "%03d" LIBXC_FUNCTIONALS_BY_NAME[substring]
                                 end
    return parse(Int, abinit_libxc_int_str)
end
function libxc_to_abinit(libxc_string::AbstractString)::Int
    psp_idx = findfirst(dict_ -> dict_["libxc"] == libxc_string, PSP8_EXCHANGE_CORRELATION)
    isnothing(psp_idx) && return libxc_to_abinit_libxc(libxc_string)
    return PSP8_EXCHANGE_CORRELATION[psp_idx]["i"]
end

save_psp(io::IO, file::Psp8File) = dict_to_psp8(io, _to_dict(file))
