"""
Build an extrapolative linear interpolator for the function `f` on a radial grid `r` where
the extrapolator returns zero outside the boundaries of the radial grid.

The various functions of the pseudopotential are symmetric around zero, so
returning zero for negative values of `r` does _not_ make physical sense.
However, all of the analysis in this package is derived for symmetric functions and does
not rely on values of `r < 0`, so these inccorect values are never used.

All of the functions decay to zero at large `r` (ideally well within the radial grid), so
extrapolating to zero _is_ valid for `r > rₙ`.
"""
function build_interpolator(f::AbstractVector{T}, r::AbstractVector{T}) where {T<:Real}
    return extrapolate(interpolate((r[firstindex(f):lastindex(f)],), f,
                                   (Gridded(Linear()),)), zero(T))
end
