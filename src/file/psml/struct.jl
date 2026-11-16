#= PSML structs: one for each element of the PSML 1.1 schema (libpsml `doc/schema/psml.rnc`),
plus the kinetic energy densities of meta-GGA pseudopotentials (PSML 1.2); each has its
element's attributes as fields (`atomic-label` is `atomic_label`) and its children after
them. Angular momenta are integers (`l="d"` is 2). Grids are where the file has them: the
data of a radial funct. is on the grid of its `radfunc`, else of its group (e.g.
`nonlocal-projectors`), else the top-level one (see `psml_grid`). =#

const PsmlAnnotation = OrderedDict{String,String}

"""
PSML `<grid>`: the radial grid, and how it was made (its annotation)
"""
struct PsmlGrid
    annotation::Union{Nothing,PsmlAnnotation}
    "Radial grid points (`grid-data`)"
    data::Vector{Float64}
    PsmlGrid(fields...) = validated(new(fields...))
end

"""
PSML `<radfunc>`: a radial function, on its own grid if it has one
"""
struct PsmlRadfunc
    grid::Union{Nothing,PsmlGrid}
    "Values (shorter than the grid if its `npts` says so)"
    data::Vector{Float64}
    PsmlRadfunc(fields...) = validated(new(fields...))
end

"""
PSML `<provenance>/<input-file>`
"""
struct PsmlInputFile
    name::String
    "The text of the input file"
    content::String
    PsmlInputFile(fields...) = validated(new(fields...))
end

"""
PSML `<provenance>`: a step of the generation of the pseudopotential
"""
struct PsmlProvenance
    record_number::Union{Nothing,Int}
    "Generation code"
    creator::String
    date::Union{Nothing,String}
    annotation::Union{Nothing,PsmlAnnotation}
    input_files::Vector{PsmlInputFile}
    PsmlProvenance(fields...) = validated(new(fields...))
end

"""
PSML `<exchange-correlation>/<libxc-info>/<functional>`
"""
struct PsmlFunctional
    "Libxc name"
    name::String
    "Libxc identifier"
    id::Int
    "`exchange`, `correlation`, `exchange-correlation` (or `XC_EXCHANGE`, ...)"
    type::Union{Nothing,String}
    weight::Union{Nothing,Float64}
    PsmlFunctional(fields...) = validated(new(fields...))
end

"""
PSML `<exchange-correlation>`
"""
struct PsmlExchangeCorrelation
    annotation::Union{Nothing,PsmlAnnotation}
    "The Libxc functionals (`libxc-info`)"
    functionals::Vector{PsmlFunctional}
    PsmlExchangeCorrelation(fields...) = validated(new(fields...))
end

"""
PSML `<shell>` of a configuration
"""
struct PsmlShell
    n::Int
    l::Int
    occupation::Float64
    occupation_up::Union{Nothing,Float64}
    occupation_down::Union{Nothing,Float64}
    PsmlShell(fields...) = validated(new(fields...))
end

"""
PSML `<valence-configuration>` or `<core-configuration>`
"""
struct PsmlConfiguration
    "`total-valence-charge` or `total-core-charge`"
    total_charge::Float64
    annotation::Union{Nothing,PsmlAnnotation}
    shells::Vector{PsmlShell}
    PsmlConfiguration(fields...) = validated(new(fields...))
end

"""
PSML `<pseudo-atom-spec>` (`<header>` before PSML 1.1)
"""
struct PsmlPseudoAtomSpec
    "Chemical symbol, or a label"
    atomic_label::String
    "Atomic number (can be fractional, for virtual atoms)"
    atomic_number::Float64
    "Pseudo-atomic (valence) charge"
    z_pseudo::Float64
    "True if there is a model core charge (non-linear core correction)"
    core_corrections::Bool
    "`no`, `scalar` or `dirac`"
    relativity::String
    "True if spin-polarized (`polarized` before PSML 1.1)"
    spin_dft::Union{Nothing,Bool}
    "True if for meta-GGA functionals (PSML 1.2)"
    meta_gga::Union{Nothing,Bool}
    "Generation scheme (e.g. `Troullier-Martins`)"
    flavor::Union{Nothing,String}
    annotation::Union{Nothing,PsmlAnnotation}
    exchange_correlation::PsmlExchangeCorrelation
    valence_configuration::PsmlConfiguration
    core_configuration::Union{Nothing,PsmlConfiguration}
    PsmlPseudoAtomSpec(fields...) = validated(new(fields...))
end

"""
PSML `<valence-charge>`: the pseudo valence charge density, with 4π prefactor
"""
struct PsmlValenceCharge
    total_charge::Float64
    is_unscreening_charge::Union{Nothing,Bool}
    rescaled_to_z_pseudo::Union{Nothing,Bool}
    annotation::Union{Nothing,PsmlAnnotation}
    radfunc::PsmlRadfunc
    PsmlValenceCharge(fields...) = validated(new(fields...))
end

"""
PSML `<pseudocore-charge>` (the model core charge density, with 4π prefactor) or
`<pseudocore-kinetic-energy-density>`
"""
struct PsmlCoreDensity
    matching_radius::Union{Nothing,Float64}
    number_of_continuous_derivatives::Union{Nothing,Int}
    annotation::Union{Nothing,PsmlAnnotation}
    radfunc::PsmlRadfunc
    PsmlCoreDensity(fields...) = validated(new(fields...))
end

"""
PSML `<valence-kinetic-energy-density>`
"""
struct PsmlValenceKineticEnergyDensity
    is_unscreening_tau::Union{Nothing,Bool}
    annotation::Union{Nothing,PsmlAnnotation}
    radfunc::PsmlRadfunc
    PsmlValenceKineticEnergyDensity(fields...) = validated(new(fields...))
end

"""
PSML `<slps>`: a semilocal potential
"""
struct PsmlSemilocalPotential
    flavor::Union{Nothing,String}
    l::Int
    j::Union{Nothing,Float64}
    n::Int
    "Core radius"
    rc::Float64
    "Reference energy"
    eref::Union{Nothing,Float64}
    radfunc::PsmlRadfunc
    PsmlSemilocalPotential(fields...) = validated(new(fields...))
end

"""
PSML `<semilocal-potentials>`: a group of semilocal potentials
"""
struct PsmlSemilocalPotentials
    "`non_relativistic`, `scalar_relativistic`, `spin_orbit`, `lj`, `up`, `down`, ..."
    set::String
    flavor::Union{Nothing,String}
    annotation::Union{Nothing,PsmlAnnotation}
    grid::Union{Nothing,PsmlGrid}
    slps::Vector{PsmlSemilocalPotential}
    PsmlSemilocalPotentials(fields...) = validated(new(fields...))
end

"""
PSML `<local-potential>`
"""
struct PsmlLocalPotential
    "How it was made (e.g. `oncv-fit`)"
    type::String
    annotation::Union{Nothing,PsmlAnnotation}
    grid::Union{Nothing,PsmlGrid}
    radfunc::PsmlRadfunc
    "The charge density whose potential is the local potential (`local-charge`)"
    local_charge::Union{Nothing,PsmlRadfunc}
    PsmlLocalPotential(fields...) = validated(new(fields...))
end

"""
PSML `<proj>`: a Kleinman-Bylander projector
"""
struct PsmlProjector
    "Kleinman-Bylander energy"
    ekb::Float64
    "Reference energy"
    eref::Union{Nothing,Float64}
    l::Int
    j::Union{Nothing,Float64}
    "Index among the projectors of the same l (and j)"
    seq::Int
    "How it was made (e.g. `oncv`, `KB`)"
    type::String
    radfunc::PsmlRadfunc
    PsmlProjector(fields...) = validated(new(fields...))
end

"""
PSML `<nonlocal-projectors>` (`<projectors>` before PSML 1.0): a group of projectors
"""
struct PsmlNonlocalProjectors
    set::String
    annotation::Union{Nothing,PsmlAnnotation}
    grid::Union{Nothing,PsmlGrid}
    projectors::Vector{PsmlProjector}
    PsmlNonlocalProjectors(fields...) = validated(new(fields...))
end

"""
PSML `<pswf>`: a pseudo-atomic wavefunction
"""
struct PsmlPseudoWaveFunction
    l::Int
    j::Union{Nothing,Float64}
    n::Int
    energy_level::Union{Nothing,Float64}
    radfunc::PsmlRadfunc
    PsmlPseudoWaveFunction(fields...) = validated(new(fields...))
end

"""
PSML `<pseudo-wave-functions>`: a group of pseudo-atomic wavefunctions
"""
struct PsmlPseudoWaveFunctions
    set::String
    annotation::Union{Nothing,PsmlAnnotation}
    grid::Union{Nothing,PsmlGrid}
    pswfs::Vector{PsmlPseudoWaveFunction}
    PsmlPseudoWaveFunctions(fields...) = validated(new(fields...))
end

"""
PSML file contents, in the layout of PSML 1.1 (and 1.2, for meta-GGA). Files of earlier
versions (0.8, 1.0) are read into it: their `header` is `pseudo_atom_spec`, the
`local-potential` and `projectors` of their `pseudopotential-operator` are
`local_potential` and `nonlocal_projectors` (its provenance is added to `provenance`), and
a text annotation is an annotation with that text under "text". They are written as
PSML 1.1. The format is described in
[García et al., Comput. Phys. Commun. 227, 51 (2018)](https://doi.org/10.1016/j.cpc.2018.02.011)
and the [PSML documentation](https://esl.cecam.org/data-standards/psml/).
"""
struct PsmlFile <: PsPFile
    "Identifier"
    identifier::String
    "PSML version of the file"
    version::String
    "`hartree`"
    energy_unit::String
    "`bohr`"
    length_unit::String
    uuid::Union{Nothing,String}
    xmlns::Union{Nothing,String}
    "The steps of the generation, in the order of the file"
    provenance::Vector{PsmlProvenance}
    pseudo_atom_spec::PsmlPseudoAtomSpec
    "The grid of the data that has no other"
    grid::Union{Nothing,PsmlGrid}
    valence_charge::PsmlValenceCharge
    pseudocore_charge::Union{Nothing,PsmlCoreDensity}
    valence_kinetic_energy_density::Union{Nothing,PsmlValenceKineticEnergyDensity}
    pseudocore_kinetic_energy_density::Union{Nothing,PsmlCoreDensity}
    semilocal_potentials::Vector{PsmlSemilocalPotentials}
    local_potential::Union{Nothing,PsmlLocalPotential}
    nonlocal_projectors::Vector{PsmlNonlocalProjectors}
    pseudo_wave_functions::Vector{PsmlPseudoWaveFunctions}
    PsmlFile(fields...) = validated(new(fields...))
end

PsmlFile(io::IO; identifier="", validation=VALIDATION[]) =
    with_validation(validation; source=identifier) do
        return _from_dict(PsmlFile, psml_to_dict(read(io, String)); identifier)
    end

function PsmlFile(path::AbstractString; identifier="", validation=VALIDATION[])
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return PsmlFile(io; identifier, validation)
    end
end

"""
The radial grid of the data of `radfunc`, in a group `group` (e.g. a
`PsmlNonlocalProjectors`, or `nothing`) of `file`: the first of their grids and that of the
file.
"""
function psml_grid(file::PsmlFile, radfunc::PsmlRadfunc, group=nothing)
    grid = psml_find_grid(file, radfunc, group)
    isnothing(grid) && error("PSML: a radial function without a grid")
    return grid.data[1:length(radfunc.data)]
end

function psml_find_grid(file::PsmlFile, radfunc::PsmlRadfunc, group)
    for x in (radfunc, group, file)
        !isnothing(x) && !isnothing(x.grid) && return x.grid
    end
    return nothing
end

#= The specification: the schema of PSML 1.1 (libpsml `doc/schema/psml.rnc`), and the
consistency of the values (e.g. of the charges). ONCVPSP writes j = 0 for l = 0, which the
schema does not have (it means j = 1/2): only a warning (`@soft_check`). =#

const PSML_SETS = ("non_relativistic", "scalar_relativistic", "spin_orbit", "lj", "up",
                   "down", "spin_average", "spin_difference", "user_extension1",
                   "user_extension2")
const PSML_FUNCTIONAL_TYPES = ("exchange", "correlation", "exchange-correlation",
                               "XC_EXCHANGE", "XC_CORRELATION", "XC_EXCHANGE_CORRELATION")
const PSML_RELATIVITY = ("no", "scalar", "dirac")

# a name without spaces (`xsd:NMTOKEN`)
psml_token(s) = !isempty(s) && !any(isspace, s)

function validate(grid::PsmlGrid)
    @check "PSML" !isempty(grid.data) && all(i -> grid.data[i] < grid.data[i+1], 1:length(grid.data)-1) "a grid is not increasing"
end

function validate(r::PsmlRadfunc)
    @check "PSML" !isempty(r.data) "a radial function without data"
    isnothing(r.grid) || @check "PSML" length(r.data) <= length(r.grid.data) "a radial function has $(length(r.data)) values on a grid of $(length(r.grid.data))"
end

validate(f::PsmlInputFile) = @check "PSML" psml_token(f.name) "input file name $(repr(f.name))"

validate(p::PsmlProvenance) =
    @check "PSML" something(p.record_number, 1) >= 1 "provenance: record-number = $(p.record_number)"

validate(f::PsmlFunctional) =
    @check "PSML" isnothing(f.type) || f.type in PSML_FUNCTIONAL_TYPES "functional $(f.name): type = $(repr(f.type)) is not one of $PSML_FUNCTIONAL_TYPES"

validate(xc::PsmlExchangeCorrelation) = @check "PSML" !isempty(xc.functionals) "no functional"

# l in s, ..., g; n in 1:9 (above l); j = l ± 1/2
function psml_check_quantum_numbers(element, l, n, j)
    @check "PSML" 0 <= l <= 4 "$element: l = $l"
    isnothing(n) || @check "PSML" 1 <= n <= 9 && n > l "$element: n = $n, l = $l"
    if !isnothing(j)
        @check "PSML" valid_j(l, j) || l == 0 && j == 0 "$element: j = $j, l = $l"
        @soft_check "PSML" !(l == 0 && j == 0) "$element: j = 0 for l = 0 (as ONCVPSP writes it), not 0.5"
    end
end

function validate(s::PsmlShell)
    psml_check_quantum_numbers("shell", s.l, s.n, nothing)
    m = 2s.l + 1  # states of each spin
    @check "PSML" 0 <= s.occupation <= 2m "shell $(s.n)$(PSML_ANGULAR_MOMENTA[s.l+1]): occupation = $(s.occupation)"
    for x in (s.occupation_up, s.occupation_down)
        isnothing(x) || @check "PSML" 0 <= x <= m "shell $(s.n)$(PSML_ANGULAR_MOMENTA[s.l+1]): spin occupation = $x"
    end
    isnothing(s.occupation_up) || isnothing(s.occupation_down) ||
        @check "PSML" isapprox(s.occupation_up + s.occupation_down, s.occupation; atol=1e-6) "shell $(s.n)$(PSML_ANGULAR_MOMENTA[s.l+1]): the spin occupations do not add to the occupation"
end

function validate(c::PsmlConfiguration)
    @check "PSML" !isempty(c.shells) "a configuration without shells"
    total = sum(s -> s.occupation, c.shells; init=0.0)
    @check "PSML" isapprox(total, c.total_charge; atol=1e-6) "a configuration of charge $(c.total_charge), but its shells hold $total"
end

function validate(spec::PsmlPseudoAtomSpec)
    @check "PSML" psml_token(spec.atomic_label) "atomic-label = $(repr(spec.atomic_label))"
    @check "PSML" spec.relativity in PSML_RELATIVITY "relativity = $(repr(spec.relativity)) is not one of $PSML_RELATIVITY"
    @check "PSML" 0 < spec.z_pseudo <= spec.atomic_number "z-pseudo = $(spec.z_pseudo), atomic-number = $(spec.atomic_number)"
    core = spec.core_configuration
    isnothing(core) ||
        @check "PSML" isapprox(spec.z_pseudo + core.total_charge, spec.atomic_number; atol=1e-6) "z-pseudo = $(spec.z_pseudo) and the core charge $(core.total_charge) do not add to atomic-number = $(spec.atomic_number)"
end

validate(c::PsmlCoreDensity) =
    @check "PSML" something(c.number_of_continuous_derivatives, 0) >= 0 "number-of-continuous-derivatives = $(c.number_of_continuous_derivatives)"

validate(p::PsmlSemilocalPotential) = psml_check_quantum_numbers("slps", p.l, p.n, p.j)
function validate(p::PsmlProjector)
    psml_check_quantum_numbers("proj", p.l, nothing, p.j)
    @check "PSML" p.seq >= 1 "proj: seq = $(p.seq)"
end
validate(w::PsmlPseudoWaveFunction) = psml_check_quantum_numbers("pswf", w.l, w.n, w.j)

# a group (e.g. of projectors): a known set, with j (only) for the set `lj`
function psml_check_group(group, name, elements)
    @check "PSML" group.set in PSML_SETS "$name: set = $(repr(group.set)) is not one of $PSML_SETS"
    @check "PSML" !isempty(elements) "$name: empty"
    @check "PSML" all(x -> isnothing(x.j) == (group.set != "lj"), elements) "$name: j is given $(group.set == "lj" ? "not for all" : "for some") of the set $(group.set)"
end
validate(g::PsmlSemilocalPotentials) = psml_check_group(g, "semilocal-potentials", g.slps)
validate(g::PsmlNonlocalProjectors) = psml_check_group(g, "nonlocal-projectors", g.projectors)
validate(g::PsmlPseudoWaveFunctions) = psml_check_group(g, "pseudo-wave-functions", g.pswfs)

# The radial functions of a file, with their group (or `nothing`), and the names
function psml_radfuncs(file::PsmlFile)
    radfuncs = Tuple{PsmlRadfunc,Any,String}[(file.valence_charge.radfunc, nothing, "valence-charge")]
    for (name, x) in ("pseudocore-charge" => file.pseudocore_charge,
                      "valence-kinetic-energy-density" => file.valence_kinetic_energy_density,
                      "pseudocore-kinetic-energy-density" => file.pseudocore_kinetic_energy_density)
        isnothing(x) || push!(radfuncs, (x.radfunc, nothing, name))
    end
    lp = file.local_potential  # a group too: it can have a grid
    if !isnothing(lp)
        push!(radfuncs, (lp.radfunc, lp, "local-potential"))
        isnothing(lp.local_charge) || push!(radfuncs, (lp.local_charge, lp, "local-charge"))
    end
    for (groups, field, name) in ((file.semilocal_potentials, :slps, "slps"),
                                  (file.nonlocal_projectors, :projectors, "proj"),
                                  (file.pseudo_wave_functions, :pswfs, "pswf"))
        for g in groups, x in getfield(g, field)
            push!(radfuncs, (x.radfunc, g, name))
        end
    end
    return radfuncs
end

function validate(file::PsmlFile)
    @check "PSML" !isnothing(tryparse(Float64, file.version)) "version = $(repr(file.version))"
    @check "PSML" file.energy_unit == "hartree" && file.length_unit == "bohr" "units $(file.energy_unit) and $(file.length_unit), not hartree and bohr"
    @check "PSML" !isempty(file.provenance) "no provenance"
    spec = file.pseudo_atom_spec
    @check "PSML" spec.core_corrections == !isnothing(file.pseudocore_charge) "core-corrections = $(spec.core_corrections), but pseudocore-charge is $(isnothing(file.pseudocore_charge) ? "missing" : "present")"
    @check "PSML" something(spec.meta_gga, false) == !isnothing(file.valence_kinetic_energy_density) "meta-gga = $(spec.meta_gga), but valence-kinetic-energy-density is $(isnothing(file.valence_kinetic_energy_density) ? "missing" : "present")"
    @check "PSML" !isempty(file.semilocal_potentials) || !isnothing(file.local_potential) "neither semilocal potentials nor a local potential"
    @check "PSML" isempty(file.nonlocal_projectors) || !isnothing(file.local_potential) "nonlocal projectors without a local potential"
    for (radfunc, group, name) in psml_radfuncs(file)
        grid = psml_find_grid(file, radfunc, group)
        @check "PSML" !isnothing(grid) "$name: no grid"
        isnothing(grid) ||
            @check "PSML" length(radfunc.data) <= length(grid.data) "$name: $(length(radfunc.data)) values on a grid of $(length(grid.data))"
    end
end

identifier(file::PsmlFile)::String = file.identifier
format(::PsmlFile)::String = "PSML"
function element(file::PsmlFile)::PeriodicTable.Element
    z = file.pseudo_atom_spec.atomic_number
    isinteger(z) && 1 <= z <= length(PeriodicTable.elements) &&
        return PeriodicTable.elements[Int(z)]
    return PeriodicTable.Element()  # a virtual atom
end
function functional(file::PsmlFile)::Vector{Symbol}
    ids = [f.id for f in file.pseudo_atom_spec.exchange_correlation.functionals]
    all(in(keys(LIBXC_FUNCTIONALS_BY_ID)), ids) || return Symbol[]  # e.g. -1: unknown
    return [Symbol(lowercase(LIBXC_FUNCTIONALS_BY_ID[id])) for id in ids]
end
valence_charge(file::PsmlFile)::Float64 = file.pseudo_atom_spec.z_pseudo
is_norm_conserving(::PsmlFile)::Bool = true
is_ultrasoft(::PsmlFile)::Bool = false
is_paw(::PsmlFile)::Bool = false
has_spin_orbit(file::PsmlFile)::Bool = file.pseudo_atom_spec.relativity == "dirac"
has_model_core_charge_density(file::PsmlFile)::Bool = !isnothing(file.pseudocore_charge)

save_psp(io::IO, file::PsmlFile) = dict_to_psml(io, _to_dict(file))
