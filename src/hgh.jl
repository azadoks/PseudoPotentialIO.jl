struct HghPsP <: AbstractPsP
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
end

function HghPsP(lines::AbstractVector{T}) where {T<:AbstractString}
    title = strip(lines[1])
    # lines[2] contains the number of electrons (and the AM channel in which they sit)
    m = match(r"^ *(([0-9]+ *)+)", lines[2])
    zion = [parse(Int, part) for part in split(m[1])]
    # lines[3] contains rloc nloc and coefficients for it
    m = match(r"^ *([-.0-9]+) +([0-9]+)( +([-.0-9]+ *)+)? *", lines[3])
    rloc = parse(Float64, m[1])
    nloc = parse(Int, m[2])
    cloc = []
    m[3] !== nothing && (cloc = [parse(Float64, part) for part in split(m[3])])
    @assert length(cloc) == nloc
    # lines[4] contains the maximal AM channel
    m = match(r"^ *([0-9]+)", lines[4])
    lmax = parse(Int, m[1]) - 1

    rp = Vector{Float64}(undef, lmax + 1)
    h = Vector{Matrix{Float64}}(undef, lmax + 1)
    cur = 5  # Current line to parse
    for l in 0:lmax
        # loop over all AM channels and extract projectors,
        # these are given in blocks like
        #
        #    0.42273813    2     5.90692831    -1.26189397
        #                                       3.25819622
        #    0.48427842    1     2.72701346
        #
        # In each such blocks, the first number is the rp, the next is the number
        # of projectors and the following numbers (including the indented continuation
        # in the next lines) is the upper triangle of the coupling matrix h for this AM.
        # Then the next block commences unindented.
        # Match the first line of a block:
        m = match(r"^ *([-.0-9]+) +([0-9]+)( +([-.0-9]+ *)+)? *", lines[cur])
        rp[l + 1] = parse(Float64, m[1])
        nproj = parse(Int, m[2])
        h[l + 1] = Matrix{Float64}(undef, nproj, nproj)
        # If there are no projectors for this AM channel nproj is zero
        # and we can just increase cur and move on to the next block.
        if nproj == 0
            cur += 1
            continue
        end
        # Else we have to parse the extra parts of the hcoeff matrix.
        # This is done here.
        hcoeff = [parse(Float64, part) for part in split(m[3])]
        for i in 1:nproj
            for j in i:nproj
                h[l + 1][j, i] = h[l + 1][i, j] = hcoeff[j - i + 1]
            end
            # Parse the next (indented) line
            cur += 1
            cur > length(lines) && break
            m = match(r"^ *(([-.0-9]+ *)+)", lines[cur])
            hcoeff = [parse(Float64, part) for part in split(m[1])]
        end
    end
    PspHgh(title, zion, rloc, nloc, cloc, lmax, rp, h)
end

HghPsP(path::AbstractString) = HghPsP(readlines(path))
HghPsP(io::IO) = HghPsP(readlines(io))

function element(psp::PspHgh)::PeriodicTable.Element
    try
        symbol = split(psp.title)[1]
        element = PeriodicTable.elements[Symbol(symbol)]
    catch
        element = "??"
    end
    return element
end
l_max(psp::PspHgh)::Int = psp.lmax
n_proj_radial(psp::PspHgh, l::Integer)::Int = length(psp.h[l+1])
n_pseudo_wfc(::PspHgh)::Int = 0
z_valence(psp::PspHgh)::Float64 = sum(psp.zion)
is_paw(::PspHgh)::Bool = false
is_ultrasoft(::PspHgh)::Bool = false
is_norm_conserving(::PspHgh)::Bool = true
is_coulomb(::PspHgh)::Bool = false
has_spin_orbit(::PspHgh)::Bool = false
has_nlcc(::PspHgh)::Bool = false
relativistic_treatment(::PspHgh)::Symbol = :scalar
format(::PspHgh)::String = "HGH"
