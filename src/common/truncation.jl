@inbounds @views function find_truncation_index(f::AbstractVector, tol)
    # eachindex -> ~3 allocations, 160 bytes
    # mean -> ~6 allocations, ~5KiB
    i_trunc = findfirst(i -> mean(abs, f[i:end]) < tol, eachindex(f))
    return isnothing(i_trunc) ? lastindex(f) : i_trunc
end

find_truncation_index(f::AbstractVector, ::Nothing) = lastindex(f)
