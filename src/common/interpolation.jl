function build_interpolator(f::AbstractVector{T}, r::AbstractVector{T}) where {T<:Real}
    return extrapolate(interpolate((r[firstindex(f):lastindex(f)],), f,
                                   (Gridded(Linear()),)), zero(T))
end
