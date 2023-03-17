"""
Build an interpolator for the function `f` on the radial grid `r` in real space and return
a function with two methods for evaluating it at radial points and arbitrary 3D vectors.

If the radial grid is linear, the function will be interpolated using cubic splines.
Otherwise, simple linear interpolation is the only method currently supported by
the interpoaltion backend.
"""
build_interpolator_real
@views @inbounds function build_interpolator_real(f::AbstractVector{T},
                                                  r::AbstractVector{T}) where {T<:Real}
    nodes = (r[firstindex(f):lastindex(f)],)
    itp = interpolate(nodes, f, (Gridded(Linear()),))
    itp_function(r) = itp(r)
    itp_function(R::AbstractVector) = itp_function(norm(R))
    return itp_function
end

@views @inbounds function build_interpolator_real(f::AbstractVector{T},
                                                  r::StepRangeLen{T}) where {T<:Real}
    itp = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))),
                r[firstindex(f):lastindex(f)])
    itp_function(r) = itp(r)
    itp_function(R::AbstractVector) = itp_function(norm(R))
    return itp_function
end
