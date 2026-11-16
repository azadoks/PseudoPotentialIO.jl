"""
Some UPF files seem to have unescaped ampersands which cause XML parsing to fail.

This could also be a problem for other XML-based formats (e.g. PSML).

In XML, `&` actually _begins_ escape sequences (e.g. `<` is `&lt;`),
so we only want to replace `&`s that are not followed by a valid escape sequence.
"""
escape_ampersand(text) = replace(text, r"&(?!(amp|apos|quot|lt|gt);)" => "&amp;")

# Numeric arrays in the text of XML elements (UPF v2, PSML)

function parse_xml_numbers(text::AbstractString)
    values = Float64[]
    for word in eachsplit(text)
        x = parse_fortran_reals(word)
        isnothing(x) && return nothing
        append!(values, x)
    end
    return values
end

function format_xml_numbers(v::AbstractVector{<:Real}; width=4)
    io = IOBuffer()
    for (i, x) in enumerate(v)
        @printf(io, " %23.16E", x)
        (i % width == 0 || i == length(v)) && print(io, '\n')
    end
    return "\n" * String(take!(io))
end

is_numeric_array(x) = x isa AbstractVector{<:Real}
