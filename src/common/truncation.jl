# Find the index where function stored in vector `f` decays to zero via a reverse
# search, assuming that the function smoothly decays to zero.
function find_truncation_index(f::AbstractVector, tol)
    i_trunc = nothing
    for (i, fi) in enumerate(Iterators.reverse(f))
        if abs(fi) > tol
            i_trunc = i
            break
        end
    end
    return isnothing(i_trunc) ? lastindex(f) : length(f) - i_trunc - 1
end

find_truncation_index(f::AbstractVector, ::Nothing) = lastindex(f)
find_truncation_index(::Nothing, ::Nothing) = nothing
