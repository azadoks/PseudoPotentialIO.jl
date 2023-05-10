"""
Build a cubic splines interpolator for the function `f` on the radial grid `r` in real
space and return a function with two methods for evaluating it at radial points and
3D vectors.
"""
build_interpolator_real
@views @inbounds function build_interpolator_real(f::AbstractVector{T},
                                                  r::AbstractVector{T}) where {T<:Real}
    itp = CubicSpline(r[firstindex(f):lastindex(f)], f)
    return itp
end
