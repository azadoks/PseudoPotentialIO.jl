"""
Build an interpolator for the function `f` on the radial grid `r` in real space and return
a function with two methods for evaluating it at radial points and arbitrary 3D vectors.
"""
function build_interpolator_real(f::AbstractVector{T}, r::AbstractVector{T}) where {T<:Real}
    itp = interpolate((r[firstindex(f):lastindex(f)],), f, (Gridded(Linear()),))
    itp_function(r) = itp(r)
    itp_function(R::AbstractVector) = itp_function(norm(R))
    return itp_function
end
