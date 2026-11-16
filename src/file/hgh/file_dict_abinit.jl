#= HGH, ABINIT formats (`pspcod` 3 and 10) <-> nested dictionaries

Based on ABINIT's reader (`psp3in` and `psp10in` in `m_psp_hgh.F90`). The layout of the
file, with the names of `HghFile` and of its `HghAbinitHeader`:
- "title": the first line
- "zatom", "zion", "pspdat": the second line
- "pspcod", "pspxc", "lmax", "lloc", "mmax", "r2well": the third line
- "rloc" and "cloc": the local part. `pspcod = 3`: the line `rloc c1 c2 c3 c4`;
  `pspcod = 10`: `rloc nloc c1 ... c_nloc`, as in the CP2K format
- "nonlocal": one entry per angular momentum, in increasing l: its radius "rp", its
  coupling coefficients "h" and, for l > 0, its spin-orbit coupling coefficients "k".
  - `pspcod = 3`: lmax + 1 entries, each the line `rp h11 h22 h33` and (l > 0) the line
    `k11 k22 k33`: "h" and "k" are the diagonals of the matrices (the rest follows from
    them); l = 3 has only `h11` and `k11`
  - `pspcod = 10`: their number on a line, then each as in the CP2K format (rp, nproj, the
    rows of the upper triangle of h), and for l > 0 the nproj rows of the upper triangle of
    k: "h" and "k" are the rows
- "trailing_text": the text after the data (as ABINIT, the angular momenta above lmax of
  `pspcod = 3` are not read)

As ABINIT does (list-directed reads), each line is read for the numbers it needs, and the
rest is ignored: the comments, and the numbers after them.
=#

"""
Read an HGH file of an ABINIT format (`pspcod` 3 or 10, its text) into nested dictionaries
(see `src/file/hgh/file_dict_abinit.jl`).
"""
function hgh_abinit_to_dict(text::AbstractString)::HghDict
    r = FortranRecords(text)
    try
        return hgh_abinit_records_to_dict(r)
    catch e
        e isa ErrorException ? error("HGH, line $(r.line): $(e.msg)") : rethrow()
    end
end

# Whether an HGH file (its text) is of an ABINIT format: its third line starts with the
# integer `pspcod` (that of the CP2K format with the real `rloc`)
function hgh_is_abinit(text::AbstractString)
    lines = split(text, '\n'; limit=4)
    length(lines) >= 3 || return false
    words = split(lines[3])
    return !isempty(words) && !isnothing(tryparse(Int, words[1]))
end

function hgh_abinit_records_to_dict(r::FortranRecords)
    d = HghDict("title" => String(strip(next_line!(r))))
    zatom, zion, pspdat = next_words!(r, 3)
    d["zatom"], d["zion"], d["pspdat"] = fortran_real(zatom), fortran_real(zion), fortran_int(pspdat)
    words = next_words!(r, 6)
    for (key, word) in zip(("pspcod", "pspxc", "lmax", "lloc", "mmax"), words)
        d[key] = fortran_int(word)
    end
    d["r2well"] = fortran_real(words[6])
    pspcod, lmax = d["pspcod"], d["lmax"]
    pspcod in (3, 10) ||
        error("pspcod = $pspcod: the HGH formats of ABINIT are pspcod 3 and 10")
    nonlocal = HghDict[]
    if pspcod == 3
        0 <= lmax <= 3 || error("lmax = $lmax is not in 0:3")
        values = hgh_reals!(r, 5)
        d["rloc"], d["cloc"] = values[1], values[2:5]
        for l in 0:lmax
            n = hgh_max_projectors(l)
            values = hgh_reals!(r, n + 1)
            channel = HghDict("rp" => values[1], "h" => values[2:end])
            l > 0 && (channel["k"] = hgh_reals!(r, n))
            push!(nonlocal, channel)
        end
    else
        words = hgh_numbers!(r, 2)
        nloc = fortran_int(words[2])
        length(words) >= nloc + 2 ||
            error("nloc = $nloc, but $(length(words) - 2) coefficients")
        d["rloc"], d["cloc"] = fortran_real(words[1]), fortran_real.(words[3:(nloc + 2)])
        for l in 0:(fortran_int(first(hgh_numbers!(r, 1))) - 1)
            words = hgh_numbers!(r, 2)
            rp, nproj = fortran_real(words[1]), fortran_int(words[2])
            length(words) >= nproj + 2 ||
                error("row 1 of h (l = $l) has $(length(words) - 2) values, but nproj = $nproj")
            h = Vector{Float64}[]
            nproj > 0 && push!(h, fortran_real.(words[3:(nproj + 2)]))
            append!(h, (hgh_reals!(r, nproj - i + 1) for i in 2:nproj))
            channel = HghDict("rp" => rp, "h" => h)
            l > 0 && (channel["k"] = Vector{Float64}[hgh_reals!(r, nproj - i + 1) for i in 1:nproj])
            push!(nonlocal, channel)
        end
    end
    d["nonlocal"] = nonlocal
    hgh_trailing_text!(d, r)
    return d
end

# The first n numbers of the next line
hgh_reals!(r, n) = fortran_real.(hgh_numbers!(r, n)[1:n])

"""
Write nested dictionaries (see `hgh_abinit_to_dict`) as an HGH file of an ABINIT format to
`io`.
"""
function dict_to_hgh_abinit(io::IO, d::AbstractDict)
    pspcod, lmax, nonlocal = d["pspcod"], d["lmax"], d["nonlocal"]
    pspcod in (3, 10) ||
        throw(ArgumentError("HGH: pspcod = $pspcod is not an HGH format (3 or 10)"))
    println(io, d["title"])
    println(io, hgh_real(d["zatom"]), hgh_real(d["zion"]), " ", lpad(string(d["pspdat"]), 6, '0'),
            "  zatom, zion, pspdat")
    println(io, join(hgh_int.((pspcod, d["pspxc"], lmax, d["lloc"], d["mmax"]))),
            hgh_real(d["r2well"]), "  pspcod, pspxc, lmax, lloc, mmax, r2well")
    if pspcod == 3
        length(d["cloc"]) == 4 ||
            throw(ArgumentError("HGH, pspcod = 3: cloc must have 4 values"))
        length(nonlocal) == lmax + 1 ||
            throw(ArgumentError("HGH, pspcod = 3: nonlocal must have lmax + 1 entries"))
        println(io, hgh_real(d["rloc"]), join(hgh_real.(d["cloc"])))
        for (l, channel) in zip(0:lmax, nonlocal)
            println(io, hgh_real(channel["rp"]), join(hgh_real.(channel["h"])))
            l > 0 && println(io, repeat(' ', 15), join(hgh_real.(channel["k"])))
        end
    else
        println(io, hgh_real(d["rloc"]), hgh_int(length(d["cloc"])), join(hgh_real.(d["cloc"])))
        println(io, hgh_int(length(nonlocal)))
        for (l, channel) in zip(0:(length(nonlocal) - 1), nonlocal)
            hgh_print_channel(io, channel)
            if l > 0
                length(channel["k"]) == length(channel["h"]) ||
                    throw(ArgumentError("HGH, pspcod = 10: k must have the rows of h"))
                hgh_print_rows(io, channel["k"], 1)
            end
        end
    end
    haskey(d, "trailing_text") && print(io, d["trailing_text"])
end
