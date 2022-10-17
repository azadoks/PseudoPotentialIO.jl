linmesh(i::Int, a::T, b::T) where {T<:Real} = a * i + b

logmesh1(i::Int, a::T, b::T) where {T<:Real} = b * exp(a * (i - 1))
logmesh1(i::Int, xmin::T, dx::T, z::T) where {T<:Real} = exp(xmin) * exp((i - 1) * dx) / z

logmesh2(i::Int, a::T, b::T) where {T<:Real} = b * (exp((i - 1) * a) - 1)
function logmesh2(i::Int, xmin::T, dx::T, z::T) where {T<:Real}
	return exp(xmin) * (exp((i - 1) * dx) - 1) / z
end

function guess_mesh_type(r::Vector{T}, rab::Vector{T}) where {T<:Real}
	nr = length(r)
	# Try linear
	a = r[2] - r[1]
	b = r[1] - a
	rguess = linmesh.(1:nr, a, b)
	if all(rguess .≈ r) & all(round.(rab, digits=4) .≈ a)
		return ("linear", a, b)
	end
	# Try log1
	a = log(r[2] / r[1])  # dx
	b = r[2] / exp(a)  # exp(xmin) / zmesh
	rguess = logmesh1.(1:nr, a, b)
	if all(rguess .≈ r) && all(rab .≈ a .* r)
		return ("log_1", a, b)
	end
	# Try log2
	b = (r[2]^2 - r[3] * r[1]) / (r[1] + r[3] - 2 * r[2])
	a = log((r[2] + b) / (r[1] + b))
	rguess = logmesh2.(1:nr, a, b)
	if all(rguess .≈ r) && all(rab .≈ a .* r .+ a * b)
		return ("log_2", a, b)
	end
	return ("unknown", NaN, NaN)
end
