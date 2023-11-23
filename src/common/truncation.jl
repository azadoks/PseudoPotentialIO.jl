# Find the index where function stored in vector `f` decays to zero via a reverse
# search, assuming that the function smoothly decays to zero.
function find_truncation_index(f::AbstractVector, tol)
    iszero(f) && return 1
    i_trunc = lastindex(f)
    for (i, fi) in enumerate(Iterators.reverse(f))
        if abs(fi) > tol
            i_trunc = i
            return lastindex(f) - i_trunc + 1
        end
    end
    return i_trunc
end

find_truncation_index(f::AbstractVector, ::Nothing) = lastindex(f)
find_truncation_index(::Nothing, _) = nothing
