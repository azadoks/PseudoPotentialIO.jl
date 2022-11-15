"""
Build an extrapolative linear interpolator for the function `f` on a radial grid `r` where
the extrapolator returns zero outside the boundaries of the radial grid.

The various functions of the pseudopotential should be reflected around the y-axis, so this
boundary condition does _not_ make physical sense for negative values of `r`. On the other
hand, all of the analysis in this package is derived for this case and does not rely on
values of `r < 0`.
All of the functions decay to zero at large `r` (ideally well within the radial grid), so
this extrapolation is motivated for `r > rₙ`.
"""
function build_interpolator(f::AbstractVector{T}, r::AbstractVector{T}) where {T<:Real}
    return extrapolate(interpolate((r[firstindex(f):lastindex(f)],), f,
                                   (Gridded(Linear()),)), zero(T))
end
