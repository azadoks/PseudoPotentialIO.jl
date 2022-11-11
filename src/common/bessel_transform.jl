@inline function bessel_transform(l::Int, r::AbstractVector{T},
                                  dr::Union{T,AbstractVector{T}}, f::AbstractVector{T},
                                  q::T)::T where {T<:Real}
    integrand = Vector{T}(undef, length(f))
    @inbounds @fastmath for i in eachindex(f)
        integrand[i] = r[i]^2 * fast_sphericalbesselj(l, q * r[i]) * f[i]
    end
    # integrand = @. r^2 * fast_sphericalbesselj(l, q * r) * f
    return 4π * trapezoid(integrand, dr)
end

@inline function bessel_transform(r::AbstractVector{T}, dr::Union{T,AbstractVector{T}},
                                  f::AbstractVector{T}, q::T)::T where {T<:Real}
    integrand = Vector{T}(undef, length(f))
    @inbounds @fastmath for i in eachindex(f)
        integrand[i] = r[i]^2 * fast_sphericalbesselj0(q * r[i]) * f[i]
    end
    # integrand = @. r^2 * fast_sphericalbesselj0(q * r) * f
    return 4π * trapezoid(integrand, dr)
end

@inline function bessel_transform(::AbstractVector{T}, ::Union{T,AbstractVector{T}},
                                  ::Nothing, ::T)::Nothing where {T<:Real}
    return nothing
end
