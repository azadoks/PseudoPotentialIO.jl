#= PSP8 <-> nested dictionaries

The layout of the file (https://docs.abinit.org/developers/psp8_info/), with the names of
`Psp8File`:
- "header": the title and the values of the header lines, by their ABINIT names (`zatom`,
  `zion`, `pspd`, `pspcod`, `pspxc`, `lmax`, `lloc`, `r2well`, `rchrg`, `fchrg`, `qchrg`,
  `extension_switch`)
- "rgrid": the radial grid (the second column of every table)
- "projectors" and "projectors_so" (spin-orbit): one block for each angular momentum with
  projectors, in the order of the file: its "l", its "ekb" and its "projectors" (a `Vector`
  of arrays, one for each column)
- "v_local": the local potential
- "nlcc": the model core charge "rhoc" and its four derivatives "d_rhoc_dr", ...
- "valence": the pseudo valence density "rhov", and the all-electron valence and core
  densities "ae_rhov" and "ae_rhoc" if the file has them (ONCVPSP writes them)
- "trailing_text": the text after the data (ONCVPSP appends its input), if any

The sizes (`mmax`, `nproj`, `nprojso`) say how much data to read; they are dropped, and
written back from the data.

Based on ABINIT's reader (`psp8in` in `m_psp8.F90`):
  - every line is read list-directed: its first values, the rest of the line is a comment
  - there is a block for each l with projectors and for `l = lloc` (the local potential),
    in increasing l, and the local potential last if `lloc > lmax`; spin-orbit blocks for
    l ≥ 1 with projectors; their first value is their l
  - the core charge is there if `fchrg > 1e-15`, the spin-orbit blocks if
    `extension_switch` is 2 or 3, the valence density if it is 1 or 3
  - the grid is read again with every table and only the last one is kept; here the first
=#

const Psp8Dict = OrderedDict{String,Any}

const PSP8_NLCC_COLUMNS = ("rhoc", "d_rhoc_dr", "d2_rhoc_dr2", "d3_rhoc_dr3", "d4_rhoc_dr4")
const PSP8_VALENCE_COLUMNS = ("rhov", "ae_rhov", "ae_rhoc")

psp8_has_nlcc(header) = header["fchrg"] > 1e-15
psp8_has_spin_orbit(header) = header["extension_switch"] in (2, 3)
psp8_has_valence(header) = header["extension_switch"] in (1, 3)

"""
Read a PSP8 file (its text) into nested dictionaries (see `src/file/psp8/file_dict.jl`).
"""
function psp8_to_dict(text::AbstractString)::Psp8Dict
    r = FortranRecords(text)
    try
        return psp8_records_to_dict(r)
    catch e
        e isa ErrorException ? error("PSP8, line $(r.line): $(e.msg)") : rethrow()
    end
end

function psp8_records_to_dict(r::FortranRecords)
    header = Psp8Dict("title" => String(next_line!(r)))
    zatom, zion, pspd = next_words!(r, 3)
    header["zatom"] = fortran_real(zatom)
    header["zion"] = fortran_real(zion)
    header["pspd"] = fortran_int(pspd)
    pspcod, pspxc, lmax, lloc, mmax, r2well = next_words!(r, 6)
    header["pspcod"] = fortran_int(pspcod)
    header["pspcod"] == 8 || error("pspcod = $(header["pspcod"]), not 8")
    header["pspxc"] = fortran_int(pspxc)
    header["lmax"] = lmax = fortran_int(lmax)
    header["lloc"] = lloc = fortran_int(lloc)
    mmax = fortran_int(mmax)
    header["r2well"] = fortran_real(r2well)
    for (key, word) in zip(("rchrg", "fchrg", "qchrg"), next_words!(r, 3))
        header[key] = fortran_real(word)
    end
    nproj = fortran_int.(next_words!(r, lmax + 1))
    lloc <= lmax && nproj[lloc + 1] > 0 &&
        error("l = lloc = $lloc has projectors (nproj = $(nproj[lloc + 1]))")
    header["extension_switch"] = fortran_int(only(next_words!(r, 1)))
    header["extension_switch"] in 0:3 ||
        error("extension_switch = $(header["extension_switch"]) is not 0, 1, 2 or 3")
    nprojso = psp8_has_spin_orbit(header) ? fortran_int.(next_words!(r, lmax)) : Int[]

    d = Psp8Dict("header" => header, "rgrid" => nothing)
    projectors = Psp8Dict[]
    v_local = nothing
    for l in 0:max(lmax, lloc)
        if l <= lmax && nproj[l + 1] > 0
            grid, block = psp8_read_projectors!(r, l, nproj[l + 1], mmax)
            push!(projectors, block)
        elseif l == lloc
            psp8_read_l!(r, l, 0)
            grid, (v_local,) = psp8_read_table!(r, mmax, 1)
        else
            continue
        end
        isnothing(d["rgrid"]) && (d["rgrid"] = grid)
    end
    d["projectors"] = projectors
    d["v_local"] = v_local
    if psp8_has_spin_orbit(header)
        d["projectors_so"] = [last(psp8_read_projectors!(r, l, nprojso[l], mmax))
                              for l in 1:lmax if nprojso[l] > 0]
    end
    if psp8_has_nlcc(header)
        _, columns = psp8_read_table!(r, mmax, 5)
        d["nlcc"] = Psp8Dict(zip(PSP8_NLCC_COLUMNS, columns))
    end
    if psp8_has_valence(header)
        # as many of the density columns as the first line has
        ncolumns = clamp(length(next_numbers!(r)) - 2, 1, 3)
        r.line -= 1
        _, columns = psp8_read_table!(r, mmax, ncolumns)
        d["valence"] = Psp8Dict(zip(PSP8_VALENCE_COLUMNS, columns))
    end
    trailing_text = remaining_text(r)
    isempty(strip(trailing_text)) || (d["trailing_text"] = trailing_text)
    return d
end

# The first line of a block: its l, then `n` values
function psp8_read_l!(r, l, n)
    words = next_words!(r, n + 1)
    fortran_int(words[1]) == l ||
        error("a block of l = $(words[1]) where l = $l is expected (the l in increasing order)")
    return fortran_real.(@view words[2:end])
end

# A block of projectors: its grid, and its dictionary
function psp8_read_projectors!(r, l, n, mmax)
    ekb = psp8_read_l!(r, l, n)
    grid, columns = psp8_read_table!(r, mmax, n)
    return grid, Psp8Dict("l" => l, "ekb" => ekb, "projectors" => columns)
end

# `mmax` lines of an index, a radius and `n` values: the grid and the `n` columns
function psp8_read_table!(r, mmax, n)
    grid = Vector{Float64}(undef, mmax)
    columns = [Vector{Float64}(undef, mmax) for _ in 1:n]
    for i in 1:mmax
        words = next_words!(r, n + 2)
        grid[i] = fortran_real(words[2])
        for (column, word) in zip(columns, @view words[3:end])
            column[i] = fortran_real(word)
        end
    end
    return grid, columns
end

"""
Write nested dictionaries (see `psp8_to_dict`) as a PSP8 file to `io`.
"""
function dict_to_psp8(io::IO, d::AbstractDict)
    header = d["header"]
    rgrid = d["rgrid"]
    lmax, lloc = header["lmax"], header["lloc"]
    blocks = psp8_checked_blocks(d, "projectors", 0:lmax)
    blocks_so = psp8_has_spin_orbit(header) ?
                psp8_checked_blocks(d, "projectors_so", 1:lmax) : Dict()
    haskey(blocks, lloc) && throw(ArgumentError("PSP8: l = lloc = $lloc has projectors"))
    for (key, present) in ("projectors_so" => psp8_has_spin_orbit(header),
                           "nlcc" => psp8_has_nlcc(header),
                           "valence" => psp8_has_valence(header))
        haskey(d, key) == present ||
            throw(ArgumentError("PSP8: `$key` does not match extension_switch and fchrg"))
    end
    nproj(blocks, l) = haskey(blocks, l) ? length(blocks[l]["projectors"]) : 0

    println(io, header["title"])
    @printf io "% 0.6f % 0.6f % 06d\n" header["zatom"] header["zion"] header["pspd"]
    @printf(io, "% 12d % 12d % 12d % 12d % 12d % 0.8f\n", header["pspcod"], header["pspxc"],
            lmax, lloc, length(rgrid), header["r2well"])
    @printf io "% 0.8f % 0.8f % 0.8f\n" header["rchrg"] header["fchrg"] header["qchrg"]
    println(io, join((@sprintf("% 12d", nproj(blocks, l)) for l in 0:lmax), ' '))
    @printf io "% 12d\n" header["extension_switch"]
    if psp8_has_spin_orbit(header)
        println(io, join((@sprintf("% 12d", nproj(blocks_so, l)) for l in 1:lmax), ' '))
    end
    for l in 0:max(lmax, lloc)
        if haskey(blocks, l)
            psp8_write_block(io, l, blocks[l]["ekb"], rgrid, blocks[l]["projectors"])
        elseif l == lloc
            psp8_write_block(io, l, nothing, rgrid, [d["v_local"]])
        end
    end
    for l in 1:lmax
        haskey(blocks_so, l) &&
            psp8_write_block(io, l, blocks_so[l]["ekb"], rgrid, blocks_so[l]["projectors"])
    end
    if psp8_has_nlcc(header)
        psp8_write_table(io, rgrid, [d["nlcc"][k] for k in PSP8_NLCC_COLUMNS])
    end
    if psp8_has_valence(header)
        psp8_write_table(io, rgrid, [d["valence"][k] for k in PSP8_VALENCE_COLUMNS
                                     if haskey(d["valence"], k)])
    end
    haskey(d, "trailing_text") && print(io, d["trailing_text"])
end

# the blocks of `d[key]` by l, checked
function psp8_checked_blocks(d, key, ls)
    blocks = Dict{Int,Any}()
    for block in d[key]
        l = block["l"]
        l in ls && !haskey(blocks, l) ||
            throw(ArgumentError("PSP8: `$key` has a block of l = $l (twice, or out of $ls)"))
        length(block["ekb"]) == length(block["projectors"]) ||
            throw(ArgumentError("PSP8: `$key` of l = $l has a different number of ekb and projectors"))
        blocks[l] = block
    end
    return blocks
end

fortran_exponential(x::Real) = replace(@sprintf("% 0.14E", x), "E" => "D")

# a block: its l and `ekb` (`nothing`: the local potential), then its table
function psp8_write_block(io, l, ekb, rgrid, columns)
    if isnothing(ekb)
        @printf io "% 19d\n" l
    else
        @printf io "% 19d %s" l repeat(" ", 21)
        foreach(e -> print(io, ' ', fortran_exponential(e)), ekb)
        println(io)
    end
    psp8_write_table(io, rgrid, columns)
end

function psp8_write_table(io, rgrid, columns)
    all(c -> length(c) == length(rgrid), columns) ||
        throw(ArgumentError("PSP8: a table that does not have the size of the grid"))
    for (i, r) in enumerate(rgrid)
        @printf io "% 19d %s" i fortran_exponential(r)
        foreach(c -> print(io, ' ', fortran_exponential(c[i])), columns)
        println(io)
    end
end
