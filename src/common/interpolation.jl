abstract type InterpolationMethod end
struct CubicSplineInterpolation <: InterpolationMethod end
struct LinearInterpolation <: InterpolationMethod end

function interpolate(m::AbstractVector{T}, y::AbstractVector{T}, ::CubicSplineInterpolation) where {T}
    x = m[eachindex(y)]
    return CubicSplines.CubicSpline(x, y)
end

function interpolate(m::AbstractVector{T}, y::AbstractVector{T}, ::LinearInterpolation) where {T}
    x = m[eachindex(y)]
    return Interpolations.interpolate((collect(x),), y, Interpolations.Gridded(Interpolations.Linear()))
end

# function interpolate(m::UniformMesh{T}, y::AbstractVector{T}, ::LinearInterpolation) where {T}
#     submesh = UniformMesh(m[eachindex(y)])
#     x = range(submesh.x1, submesh.xn, submesh.n)
#     return Interpolations.linear_interpolation(x, y)
# end

function interpolate_evaluate(m::AbstractVector{T}, y::AbstractVector{T}, xnew::AbstractVector{T}, interpolation_method::InterpolationMethod) where{T}
    itp = interpolate(m, y, interpolation_method)
    return itp.(xnew)
end

function interpolate(psp::NumericPsP, quantity::AbstractPsPQuantity, args...; interpolation_method=CubicSplineInterpolation())
    !has_quantity(psp, quantity) && return Nothing
    y = get_quantity(psp, quantity, args...)
    return interpolate(psp.r, y, interpolation_method)
end

function interpolate_evaluate(psp::NumericPsP, mesh_new::AbstractVector, quantity::AbstractPsPQuantity, args...; interpolation_method=CubicSplineInterpolation())
    !has_quantity(psp, quantity) && return Nothing
    y = get_quantity(psp, quantity, args...)
    return interpolate_evaluate(psp.r, y, mesh_new, interpolation_method)
end

# function CubicSplines.CubicSpline(psp::NumericPsP, quantity::AbstractPsPQuantity, args...)
#     !has_quantity(psp, quantity) && return nothing
#     f = get_quantity(psp, quantity, args...)
#     return CubicSplines.CubicSpline(psp.r[eachindex(f)], f)
# end

# function cubic_spline_interpolate(psp::NumericPsP, mesh::RadialMesh, quantity::AbstractPsPQuantity, args...)
#     itp = CubicSpline(psp, quantity, args...)
#     isnothing(itp) && return nothing
#     rmax = cutoff_radius(psp, quantity, args...)
#     imax = index_leq(mesh, rmax)
#     return itp.(mesh[1:imax])
# end

function interpolate_evaluate(psp::NumericPsP, mesh::RadialMesh, quantity::AbstractProjector; interpolation_method=CubicSplineInterpolation())
    if has_quantity(psp, quantity)
        f = map(angular_momenta(psp)) do l
            map(1:n_radials(psp, quantity, l)) do n
                return interpolate_evaluate(psp, mesh, quantity, l, n; interpolation_method)
            end
        end
        f = OffsetArray(f, angular_momenta(psp))
    else
        f = nothing
    end

    return f
end

function interpolate_evaluate(psp::NormConservingPsP{T},
                                  mesh::RadialMesh{T}; interpolation_method=CubicSplineInterpolation())::NormConservingPsP{T} where {T}
    Vloc = interpolate_evaluate(psp, mesh, NumericLocalPotential(); interpolation_method)

    β = interpolate_evaluate(psp, mesh, NumericProjector(); interpolation_method)
    χ = interpolate_evaluate(psp, mesh, NumericState(); interpolation_method)
    ρval = interpolate_evaluate(psp, mesh, ValenceChargeDensity(); interpolation_method)
    ρcore = interpolate_evaluate(psp, mesh, CoreChargeDensity(); interpolation_method)

    identifier = "$(psp.identifier) on $(mesh)"
    checksum = nothing

    return NormConservingPsP{T}(identifier, checksum, psp.Zatom, psp.Zval, psp.lmax, mesh,
                                Vloc, β, psp.D, χ, ρcore, ρval)
end

function interpolate_evaluate(psp::UltrasoftPsP{T},
                                  mesh::RadialMesh{T}; interpolation_method=CubicSplineInterpolation())::UltrasoftPsP{T} where {T}
    Vloc = interpolate_evaluate(psp, mesh, NumericLocalPotential(); interpolation_method)

    β = interpolate_evaluate(psp, mesh, NumericProjector(); interpolation_method)
    χ = interpolate_evaluate(psp, mesh, NumericState(); interpolation_method)
    ρval = interpolate_evaluate(psp, mesh, ValenceChargeDensity(); interpolation_method)
    ρcore = interpolate_evaluate(psp, mesh, CoreChargeDensity(); interpolation_method)

    Q = map(angular_momenta(psp)) do l
        Ql = similar(psp.Q[l])
        for n in axes(psp.Q[l], 1)
            for m in axes(psp.Q[l], 2)
                Ql[n, m] = interpolate_evaluate(psp, mesh, AugmentationFunction(), l, n,
                                                    m; interpolation_method)
            end
        end
        return Ql
    end
    Q = OffsetVector(Q, angular_momenta(psp))

    identifier = "$(psp.identifier) on $(mesh)"
    checksum = nothing

    return UltrasoftPsP{T}(identifier, checksum, psp.Zatom, psp.Zval, psp.lmax, mesh, Vloc,
                           β, psp.D, χ, Q, psp.q, ρcore, ρval)
end
