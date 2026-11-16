#= HGH <-> nested dictionaries

The layout of the file, with the names of `HghFile`:
- "title": the first line
- "zion": the number of valence electrons of each angular momentum (second line)
- "rloc" and "cloc": the radius and the coefficients of the local part (third line)
- "nonlocal": one entry per angular momentum (their number is the fourth line), in
  increasing l: its radius "rp" and the rows of the upper triangle of its coupling matrix
  "h" (a `Vector` of arrays, one for each line: row i has nproj - i + 1 values)
- "trailing_text": the text after the data

The sizes (`nloc`, the number of angular momenta, of projectors) say how much data to read
and are checked against the values on their lines; they are dropped, and written back from
the data. A line is read up to its first word that is not a number: the rest is a comment.
=#

const HghDict = OrderedDict{String,Any}

"""
Read an HGH file (its text) into nested dictionaries (see `src/file/hgh/file_dict.jl`).
"""
function hgh_to_dict(text::AbstractString)::HghDict
    r = FortranRecords(text)
    try
        return hgh_records_to_dict(r)
    catch e
        e isa ErrorException ? error("HGH, line $(r.line): $(e.msg)") : rethrow()
    end
end

function hgh_records_to_dict(r::FortranRecords)
    d = HghDict("title" => String(strip(next_line!(r))))
    d["zion"] = fortran_int.(hgh_numbers!(r, 1))
    words = hgh_numbers!(r, 2)
    d["rloc"] = fortran_real(words[1])
    nloc = fortran_int(words[2])
    length(words) == nloc + 2 ||
        error("nloc = $nloc, but $(length(words) - 2) coefficients")
    d["cloc"] = fortran_real.(words[3:end])
    nonlocal = HghDict[]
    for l in 0:(fortran_int(first(hgh_numbers!(r, 1))) - 1)
        words = hgh_numbers!(r, 2)
        rp, nproj = fortran_real(words[1]), fortran_int(words[2])
        values = @view words[3:end]
        h = Vector{Float64}[]
        for i in 1:nproj
            i > 1 && (values = next_numbers!(r))
            length(values) == nproj - i + 1 ||
                error("row $i of h (l = $l) has $(length(values)) values, but nproj = $nproj")
            push!(h, fortran_real.(values))
        end
        nproj == 0 && !isempty(values) && error("h (l = $l) has values, but nproj = 0")
        push!(nonlocal, HghDict("rp" => rp, "h" => h))
    end
    d["nonlocal"] = nonlocal
    trailing_text = remaining_text(r)
    isempty(strip(trailing_text)) || (d["trailing_text"] = trailing_text)
    return d
end

function hgh_numbers!(r, n)
    words = next_numbers!(r)
    length(words) >= n || error("$(length(words)) numbers, but at least $n are needed")
    return words
end

"""
Write nested dictionaries (see `hgh_to_dict`) as an HGH file to `io`.
"""
function dict_to_hgh(io::IO, d::AbstractDict)
    println(io, d["title"])
    println(io, join(hgh_int.(d["zion"])))
    println(io, hgh_real(d["rloc"]), hgh_int(length(d["cloc"])), join(hgh_real.(d["cloc"])))
    println(io, hgh_int(length(d["nonlocal"])))
    for channel in d["nonlocal"]
        h = channel["h"]
        nproj = length(h)
        all(length(h[i]) == nproj - i + 1 for i in 1:nproj) ||
            throw(ArgumentError("HGH: the rows of h must have nproj, nproj - 1, ... values"))
        print(io, hgh_real(channel["rp"]), hgh_int(nproj))
        for (i, row) in enumerate(h)  # under their columns
            i > 1 && print(io, repeat(' ', 20 + 15 * (i - 1)))
            println(io, join(hgh_real.(row)))
        end
        nproj == 0 && println(io)
    end
    haskey(d, "trailing_text") && print(io, d["trailing_text"])
end

hgh_real(x) = " " * lpad(string(Float64(x)), 14)  # shortest text that reads back the same
hgh_int(n) = " " * lpad(string(n), 4)
