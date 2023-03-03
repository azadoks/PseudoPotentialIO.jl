@inbounds @views function find_truncation_index(f::AbstractVector, tol)
    i_trunc = findfirst(i -> mean(abs, f[i:end]) < tol, eachindex(f))
    return isnothing(i_trunc) ? lastindex(f) : i_trunc
end

find_truncation_index(f::AbstractVector, ::Nothing) = lastindex(f)

find_truncation_index(::Nothing, ::Nothing) = nothing
