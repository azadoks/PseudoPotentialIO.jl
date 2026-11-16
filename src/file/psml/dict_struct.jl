#= PSML dictionaries (`file_dict.jl`) <-> `PsmlFile`: `_from_dict` builds the structs (from
the layout of PSML 1.1, or of earlier versions: see `PsmlFile`), `_to_dict` is its inverse
and gives the layout of PSML 1.1. A struct is its element's attributes field for field
(`psml_struct`, `psml_fields`), then its children. =#

const PSML_ANGULAR_MOMENTA = "spdfghi"

function psml_l(l)
    i = l isa AbstractString && length(l) == 1 ? findfirst(only(l), PSML_ANGULAR_MOMENTA) : nothing
    isnothing(i) && error("PSML: l = $(repr(l)) is not one of $(join(PSML_ANGULAR_MOMENTA, ", "))")
    return i - 1
end

# the attribute of a field: `atomic_label` is `atomic-label`, but the units and
# `energy_level` are spelled with underscores
psml_key(f::Symbol) = f in (:energy_unit, :length_unit, :energy_level) ? string(f) :
                      replace(string(f), '_' => '-')

# Entry `key` of `d` as a `T`: a missing entry is `nothing` if `T` allows it, an error if not.
function psml_get(::Type{T}, d, key::AbstractString) where {T}
    S = Base.nonnothingtype(T)
    x = get(d, key, nothing)
    if isnothing(x)
        Nothing <: T && return nothing
        error("PSML: `$key` is missing")
    end
    try
        return S === String ? string(x) : convert(S, x)
    catch
        error("PSML: cannot read `$key` = $(repr(x)) as $S")
    end
end

# an annotation: its attributes (before PSML 1.0, it could be a text attribute)
psml_annotation(::Nothing) = nothing
psml_annotation(a::AbstractString) = PsmlAnnotation("text" => a)
psml_annotation(a::AbstractDict) = PsmlAnnotation(k => string(v) for (k, v) in a)

psml_from(::Type{PsmlGrid}, d::AbstractDict) =
    PsmlGrid(psml_annotation(get(d, "annotation", nothing)),
             psml_get(Vector{Float64}, d, "grid-data"))

psml_from(::Type{PsmlRadfunc}, d::AbstractDict) =
    PsmlRadfunc(psml_optional(PsmlGrid, d, "grid"), psml_get(Vector{Float64}, d, "data"))

# an optional element `key` of `d` as a `T`
psml_optional(::Type{T}, d, key) where {T} = haskey(d, key) ? psml_from(T, d[key]) : nothing
psml_from(::Type{T}, d::AbstractDict) where {T} = psml_struct(T, d)

psml_radfunc(d) = psml_from(PsmlRadfunc, psml_get(AbstractDict, d, "radfunc"))

# A field of a struct from its element `d`
function psml_field(::Type{T}, d, f::Symbol) where {T}
    f === :annotation && return psml_annotation(get(d, "annotation", nothing))
    f === :grid && return psml_optional(PsmlGrid, d, "grid")
    f === :radfunc && return psml_radfunc(d)
    f === :l && return psml_l(psml_get(String, d, "l"))
    return psml_get(fieldtype(T, f), d, psml_key(f))
end

# A struct from its element `d`, field by field; `fields` gives other values.
psml_struct(::Type{T}, d; fields...) where {T} =
    T((haskey(fields, f) ? fields[f] : psml_field(T, d, f) for f in fieldnames(T))...)

# the elements `key` of `d`, as a vector
psml_list(d, key) = (x = get(d, key, Any[]); x isa AbstractVector ? x : [x])

function psml_configuration(d, total)
    shells = [psml_struct(PsmlShell, s) for s in psml_list(d, "shell")]
    return psml_struct(PsmlConfiguration, d; total_charge=psml_get(Float64, d, total), shells)
end

function _from_dict(::Type{PsmlPseudoAtomSpec}, d::AbstractDict)
    xc = psml_get(AbstractDict, d, "exchange-correlation")
    libxc = psml_get(AbstractDict, xc, "libxc-info")
    functionals = [psml_struct(PsmlFunctional, f) for f in psml_list(libxc, "functional")]
    n = psml_get(Int, libxc, "number-of-functionals")
    @check "PSML" length(functionals) == n "number-of-functionals = $n, but $(length(functionals)) <functional>"
    core = get(d, "core-configuration", nothing)
    return psml_struct(PsmlPseudoAtomSpec, d;
                       core_corrections=something(get(d, "core-corrections", nothing), false),
                       spin_dft=get(d, "spin-dft", get(d, "polarized", nothing)),  # < 1.1
                       exchange_correlation=PsmlExchangeCorrelation(
                           psml_annotation(get(xc, "annotation", nothing)), functionals),
                       valence_configuration=psml_configuration(
                           psml_get(AbstractDict, d, "valence-configuration"),
                           "total-valence-charge"),
                       core_configuration=isnothing(core) ? nothing :
                                          psml_configuration(core, "total-core-charge"))
end

function _from_dict(::Type{PsmlFile}, d::AbstractDict; identifier="")
    provenance = [psml_struct(PsmlProvenance, p;
                              input_files=[psml_struct(PsmlInputFile, f)
                                           for f in psml_list(p, "input-file")])
                  for p in psml_list(d, "provenance")]
    # before PSML 1.0: the local potential and projectors are in an operator
    operator = get(d, "pseudopotential-operator", d)
    if operator !== d
        append!(provenance, psml_struct(PsmlProvenance, p; input_files=PsmlInputFile[])
                            for p in psml_list(operator, "provenance"))
    end
    spec = get(d, "pseudo-atom-spec", get(d, "header", nothing))  # header: < 1.1
    isnothing(spec) && error("PSML: <pseudo-atom-spec> is missing")
    local_potential = get(operator, "local-potential", nothing)
    if !isnothing(local_potential)
        local_charge = get(local_potential, "local-charge", nothing)
        local_potential = psml_struct(PsmlLocalPotential, local_potential;
                                      local_charge=isnothing(local_charge) ? nothing :
                                          psml_radfunc(local_charge))
    end
    nonlocal = vcat(psml_list(operator, "nonlocal-projectors"),
                    psml_list(operator, "projectors"))  # projectors: < 1.0
    return psml_struct(PsmlFile, d;
                       identifier, provenance,
                       pseudo_atom_spec=_from_dict(PsmlPseudoAtomSpec, spec),
                       valence_charge=psml_struct(PsmlValenceCharge,
                                                  psml_get(AbstractDict, d, "valence-charge")),
                       pseudocore_charge=psml_optional(PsmlCoreDensity, d, "pseudocore-charge"),
                       valence_kinetic_energy_density=psml_optional(
                           PsmlValenceKineticEnergyDensity, d, "valence-kinetic-energy-density"),
                       pseudocore_kinetic_energy_density=psml_optional(
                           PsmlCoreDensity, d, "pseudocore-kinetic-energy-density"),
                       semilocal_potentials=[psml_struct(PsmlSemilocalPotentials, g;
                                                         slps=[psml_struct(PsmlSemilocalPotential, p)
                                                               for p in psml_list(g, "slps")])
                                             for g in psml_list(d, "semilocal-potentials")],
                       local_potential,
                       nonlocal_projectors=[psml_struct(PsmlNonlocalProjectors, g;
                                                        projectors=[psml_struct(PsmlProjector, p)
                                                                    for p in psml_list(g, "proj")])
                                            for g in nonlocal],
                       pseudo_wave_functions=[psml_struct(PsmlPseudoWaveFunctions, g;
                                                          pswfs=[psml_struct(PsmlPseudoWaveFunction, p)
                                                                 for p in psml_list(g, "pswf")])
                                              for g in psml_list(d, "pseudo-wave-functions")])
end

# A value of a struct as a dictionary entry
psml_value(::Symbol, x) = x
psml_value(f::Symbol, l::Int) = f === :l ? string(PSML_ANGULAR_MOMENTA[l + 1]) : l
function psml_value(::Symbol, grid::PsmlGrid)
    d = PsmlDict()
    isnothing(grid.annotation) || (d["annotation"] = grid.annotation)
    d["grid-data"] = grid.data
    return d
end
function psml_value(::Symbol, radfunc::PsmlRadfunc)
    d = PsmlDict()
    isnothing(radfunc.grid) || (d["grid"] = psml_value(:grid, radfunc.grid))
    d["data"] = radfunc.data
    return d
end

# The element of a struct: its fields that are not `nothing`, but those of `skip`
psml_fields(x, skip...) =
    PsmlDict(psml_key(f) => psml_value(f, getfield(x, f)) for f in fieldnames(typeof(x))
             if !(f in skip) && !isnothing(getfield(x, f)))

function psml_configuration_dict(c::PsmlConfiguration, total)
    d = PsmlDict(total => c.total_charge)
    isnothing(c.annotation) || (d["annotation"] = c.annotation)
    d["shell"] = psml_fields.(c.shells)
    return d
end

function _to_dict(spec::PsmlPseudoAtomSpec)
    d = psml_fields(spec, :exchange_correlation, :valence_configuration,
                        :core_configuration)
    xc = spec.exchange_correlation
    d["exchange-correlation"] = xc_dict = PsmlDict()
    isnothing(xc.annotation) || (xc_dict["annotation"] = xc.annotation)
    xc_dict["libxc-info"] = PsmlDict("number-of-functionals" => length(xc.functionals),
                                     "functional" => psml_fields.(xc.functionals))
    d["valence-configuration"] = psml_configuration_dict(spec.valence_configuration,
                                                         "total-valence-charge")
    isnothing(spec.core_configuration) ||
        (d["core-configuration"] = psml_configuration_dict(spec.core_configuration,
                                                           "total-core-charge"))
    return d
end

# a group (e.g. of projectors), with its elements (if any) under `key`
function psml_group(group, field, key)
    d = psml_fields(group, field)
    elements = getfield(group, field)
    isempty(elements) || (d[key] = psml_fields.(elements))
    return d
end

function _to_dict(file::PsmlFile)
    # earlier versions are written in the layout of 1.1
    version = VersionNumber(file.version) < v"1.1" ? "1.1" : file.version
    d = PsmlDict("version" => version)
    for f in (:energy_unit, :length_unit, :uuid, :xmlns)
        isnothing(getfield(file, f)) || (d[psml_key(f)] = getfield(file, f))
    end
    d["provenance"] = [psml_group(p, :input_files, "input-file") for p in file.provenance]
    d["pseudo-atom-spec"] = _to_dict(file.pseudo_atom_spec)
    isnothing(file.grid) || (d["grid"] = psml_value(:grid, file.grid))
    d["valence-charge"] = psml_fields(file.valence_charge)
    for f in (:pseudocore_charge, :valence_kinetic_energy_density,
              :pseudocore_kinetic_energy_density)
        isnothing(getfield(file, f)) || (d[psml_key(f)] = psml_fields(getfield(file, f)))
    end
    d["semilocal-potentials"] = [psml_group(g, :slps, "slps") for g in file.semilocal_potentials]
    if !isnothing(file.local_potential)
        d["local-potential"] = local_potential = psml_fields(file.local_potential, :local_charge)
        lc = file.local_potential.local_charge
        isnothing(lc) ||
            (local_potential["local-charge"] = PsmlDict("radfunc" => psml_value(:radfunc, lc)))
    end
    d["nonlocal-projectors"] = [psml_group(g, :projectors, "proj")
                                for g in file.nonlocal_projectors]
    d["pseudo-wave-functions"] = [psml_group(g, :pswfs, "pswf")
                                  for g in file.pseudo_wave_functions]
    for key in ("semilocal-potentials", "nonlocal-projectors", "pseudo-wave-functions")
        isempty(d[key]) && delete!(d, key)
    end
    return d
end
