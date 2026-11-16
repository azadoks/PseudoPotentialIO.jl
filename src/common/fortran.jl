# Fortran list-directed input (PSP8, HGH, and the numbers of UPF v2 and PSML)

# A Fortran list-directed real: `D` exponents, an exponent without its letter
# (`2.0-319`), `r*value` repeats. `nothing` if `word` is not one.
function parse_fortran_reals(word::AbstractString)
    count, value = occursin('*', word) ? split(word, '*'; limit=2) : ("1", word)
    value = replace(value, r"[dD]" => "E")
    value = replace(value, r"(?<=[0-9.])([+-][0-9]+)$" => s"E\1")
    n = tryparse(Int, count)
    x = tryparse(Float64, value)
    return isnothing(n) || isnothing(x) ? nothing : fill(x, n)
end

function fortran_real(word::AbstractString)
    x = parse_fortran_reals(word)
    (isnothing(x) || length(x) != 1) && error("Not a real number: $(repr(word))")
    return only(x)
end

function fortran_int(word::AbstractString)
    x = tryparse(Int, word)
    isnothing(x) && error("Not an integer: $(repr(word))")
    return x
end

# The lines of a text, read one at a time like Fortran records
mutable struct FortranRecords
    lines::Vector{SubString{String}}
    line::Int  # the last line read
end
FortranRecords(text::AbstractString) = FortranRecords(split(String(text), '\n'), 0)

function next_line!(r::FortranRecords)
    r.line < length(r.lines) || error("Line $(r.line + 1): the file ends")
    r.line += 1
    return r.lines[r.line]
end

# The first `n` words of the next line: a list-directed read ignores the rest of the line
# (the comments of PSP8 headers, `zatom,zion,pspd`)
function next_words!(r::FortranRecords, n::Integer)
    words = split(next_line!(r))
    length(words) >= n ||
        error("Line $(r.line): $(length(words)) values, but $n are needed")
    return @view words[1:n]
end

# The words of the next line up to the first one that is not a number (a comment)
function next_numbers!(r::FortranRecords)
    words = split(next_line!(r))
    n = something(findfirst(w -> isnothing(parse_fortran_reals(w)), words), length(words) + 1)
    return @view words[1:(n - 1)]
end

# The lines that have not been read
remaining_text(r::FortranRecords) = join(@view(r.lines[(r.line + 1):end]), '\n')
