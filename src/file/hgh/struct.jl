"""
Hartwigsen-Goedecker-Hutter pseudopotential file contents.
"""
struct HghFile <: PsPFile
    "Identifier"
    identifier::String
    "Functional guessed from the source file path if available"
    functional::String
    "Description"
    title::String
    "Pseudo-atomic (valence) charge"
    zion::Vector{Int}
    "Cutoff radius for the local part of the pseudopotential"
    rloc::Float64
    "Number of coefficients defining the local part of the pseudopotential"
    nloc::Int
    "Coefficients of the local part of the pseudopotential"
    cloc::Vector{Float64}
    "Maximum angular momentum"
    lmax::Int
    "Non-local projector cutoff radius for each angular momentum"
    rp::Vector{Float64}
    "Kleinman-Bylander energies"
    h::Vector{Matrix{Float64}}
    HghFile(fields...) = validated(new(fields...))
end

function validate(file::HghFile)
    @check "HGH" all(>=(0), file.zion) && sum(file.zion) > 0 "zion = $(file.zion)"
    @check "HGH" file.rloc > 0 "rloc = $(file.rloc)"
    @check "HGH" file.nloc == length(file.cloc) "nloc = $(file.nloc), but $(length(file.cloc)) coefficients"
    @check "HGH" 0 <= file.nloc <= 4 "nloc = $(file.nloc) is not in 0:4"
    @check "HGH" length(file.rp) == length(file.h) == file.lmax + 1 "lmax = $(file.lmax), but $(length(file.rp)) radii and $(length(file.h)) coupling matrices"
    for (l, (rp, h)) in enumerate(zip(file.rp, file.h))
        @check "HGH" rp > 0 "rp = $rp for l = $(l - 1)"
        @check "HGH" size(h, 1) == size(h, 2) <= 3 && h == permutedims(h) "the coupling matrix of l = $(l - 1) is not symmetric of size ≤ 3"
    end
end

HghFile(io::IO; identifier="", validation=VALIDATION[]) =
    with_validation(validation; source=identifier) do
        return _from_dict(HghFile, hgh_to_dict(read(io, String)); identifier)
    end

function HghFile(path::AbstractString; identifier="", validation=VALIDATION[])
    identifier = isempty(identifier) ? splitpath(path)[end] : identifier
    open(path, "r") do io
        return HghFile(io; identifier, validation)
    end
end

identifier(psp::HghFile)::String = psp.identifier
format(::HghFile)::String = "HGH"
function element(psp::HghFile)::PeriodicTable.Element
    # Guess the element from the title line, if possible
    title = split(psp.title)
    isempty(title) && return PeriodicTable.Element()
    symbol = title[1]
    return haskey(PeriodicTable.elements, Symbol(symbol)) ? PeriodicTable.elements[Symbol(symbol)] : PeriodicTable.Element()
end
functional(psp::HghFile)::Vector{Symbol} = map(name -> Symbol(lowercase(name)), split(psp.functional))
valence_charge(psp::HghFile)::Float64 = sum(psp.zion, init=0)
is_norm_conserving(::HghFile)::Bool = true
is_ultrasoft(::HghFile)::Bool = false
is_paw(::HghFile)::Bool = false
has_spin_orbit(::HghFile)::Bool = false
has_model_core_charge_density(::HghFile)::Bool = false

save_psp(io::IO, file::HghFile) = dict_to_hgh(io, _to_dict(file))
