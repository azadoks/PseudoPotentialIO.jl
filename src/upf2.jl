using EzXML

export parse_upf2

parse_bool(s::AbstractString) = occursin("T", uppercase(s)) ? true : false
parse_bool(s::Char) = uppercase(s) == 'T' ? true : false

function get_attr(::Type{T}, node::EzXML.Node, key; default=zero(T)) where {T<:Number}
	return haskey(node, key) ? parse(T, strip(node[key])) : default
end

function get_attr(::Type{T}, node::EzXML.Node, key; default=T("")) where {T<:AbstractString}
	return haskey(node, key) ? T(strip(node[key])) : default
end

function get_attr(::Type{Bool}, node::EzXML.Node, key; default=false)
	return haskey(node, key) ? parse_bool(strip(node[key])) : default
end

function get_content(::Type{T}, node::EzXML.Node, dims...) where {T<:Number}
	text = split(strip(node.content))
	@assert length(text) == prod(dims)
	value = Array{T}(undef, dims...)
	for i in eachindex(value)
		value[i] = parse(T, text[i])
	end
	return value
end

function get_content(::Type{T}, node::EzXML.Node, key, dims...) where {T<:Number}
	target_node = findfirst(node, key)
	content = ismissing(target_node) ? nothing : get_content(T, target_node, dims...)
	return content
end

"""
	parse_header_upf2!(doc_root::EzXML.Node, upf::Dict)

Parse header (`PP_HEADER`) data, storing it in `upf["header"]::Dict` with the following
contents:
- `format_version::Int`: always 2
- `generated::String`: code that generated the pseudo
- `author::String`
- `date::String`: generation date (arbitrary format)
- `comment::String`
- `element::String`: elemental symbol
- `pseudo_type::String`: `NC` (norm-conserving), `US` (ultrasoft), `PAW`
(plane-augmented wave) `SL` (semi-local), `1/r` (Coulomb)
- `relativistic::String`: `scalar`, `full`, or `nonrelativistic`
- `is_ultrasoft::Bool`
- `is_paw::Bool`
- `is_coulomb::Bool`: fake Coulomb potential for all-electron calculations
- `has_so::Bool`: has spin-orbit coupling
- `has_wfc::Bool`: has pseudo-atomic wavefunctions
- `has_gipaw::Bool`: has GIPAW reconstruction data
- `paw_as_gipaw::Bool`: ???
- `core_correction::Bool`: non-linear core-correction
- `functional::String`
- `z_valence::Float64`: pseudo-ion charge
- `total_psenergy::Float64`: total energy of the pseudo-ion
- `ecutwfc::Float64`: recommended plane-wave energy cutoff
- `ecutrho::Float64`: recommended charge-density energy cutoff
- `l_max::Int`: maximum angular momentum channel of the Kleinman-Bylander projectors
- `l_max_rho::Int`
- `l_local::Int`: angular momentum channel of the local potential (-1 if none)
- `mesh_size::Int`: number of points in the radial mesh
- `number_of_wfc::Int`: number of pseudo-atomic wavefunctions
- `number_of_proj::Int`: number of Kleinman-Bylander projectors
"""
function parse_header_upf2!(doc_root::EzXML.Node, upf::Dict)
	header = Dict()
	# version = get_attr(String, doc_root, "version")
	node = findfirst("PP_HEADER", doc_root)
	header["format_version"] = 2
	header["generated"] = get_attr(String, node, "generated")
	header["author"] = get_attr(String, node, "author")
	header["date"] = get_attr(String, node, "date")
	header["comment"] = get_attr(String, node, "comment")
	header["element"] = get_attr(String, node, "element")
	header["pseudo_type"] = get_attr(String, node, "pseudo_type")
	header["relativistic"] = get_attr(String, node, "relativistic")
	header["is_ultrasoft"] = get_attr(Bool, node, "is_ultrasoft")
	header["is_paw"] = get_attr(Bool, node, "is_paw")
	header["is_coulomb"] = get_attr(Bool, node, "is_coulomb")
	header["has_so"] = get_attr(Bool, node, "has_so")
	header["has_wfc"] = get_attr(Bool, node, "has_wfc")
	header["has_gipaw"] = get_attr(Bool, node, "has_gipaw")
	header["paw_as_gipaw"] = get_attr(Bool, node, "paw_as_gipaw")
	header["core_correction"] = get_attr(Bool, node, "core_correction")
	header["functional"] = get_attr(String, node, "functional")
	header["z_valence"] = get_attr(Float64, node, "z_valence")
	header["total_psenergy"] = get_attr(Float64, node, "total_psenergy")
	header["wfc_cutoff"] = get_attr(Float64, node, "wfc_cutoff")
	header["rho_cutoff"] = get_attr(Float64, node, "rho_cutoff")
	header["l_max"] = get_attr(Int, node, "l_max")
	header["l_max_rho"] = get_attr(Int, node, "l_max_rho")
	header["l_local"] = get_attr(Int, node, "l_local")
	header["mesh_size"] = get_attr(Int, node, "mesh_size")
	header["number_of_wfc"] = get_attr(Int, node, "number_of_wfc")
	header["number_of_proj"] = get_attr(Int, node, "number_of_proj")

	return upf["header"] = header
end

"""
	parse_radial_grid_upf2!(io::IO, upf::Dict)

Parse radial grid data (`<PP_R>`) and integration factors (`<PP_RAB>`) from the `<PP_MESH>`
block, storing them in `upf["radial_grid"]` and `upf["radial_grid_derivative"]`
respectively. Also, parse the attributes of `<PP_MESH>` and store them in
`upf["radial_grid_parameters"]`:
- `dx::Float64'
- `mesh::Int`: number of mesh points
- `xmin::Float64`
- `rmax::Float64`
- `zmesh::Float64`
If these mesh parameters are present, the radial grid is one of the following:
"log_1"
``e^{x_\\text{min}} e^{(i - 1)dx} / Z_\\text{mesh}``
"log_2"
``e^{x_\\text{min}} (e^{(i - 1)dx} - 1) / Z_\\text{mesh}``
Otherwised, it is likely a linear mesh.
The type of mesh is stored in `upf["radial_grid_parameters"]["mesh_type"]`.

The radial grid derivative is the factor for discrete integration:
``\\int f(r) dr = \\sum_{i=1}^{N} f(i) r_{ab}(i)``  
"""
function parse_radial_grid_upf2!(doc_root::EzXML.Node, upf::Dict)
	node = findfirst("PP_MESH/PP_R", doc_root)
	upf["radial_grid"] = parse.(Float64, split(strip(node.content)))  # Bohr

	node = findfirst("PP_MESH/PP_RAB", doc_root)
	upf["radial_grid_derivative"] = parse.(Float64, split(strip(node.content)))

	node = findfirst("PP_MESH", doc_root)
	dx = get_attr(Float64, node, "dx")
	mesh = get_attr(Int, node, "mesh")
	xmin = get_attr(Float64, node, "xmin")
	rmax = get_attr(Float64, node, "rmax")
	zmesh = get_attr(Float64, node, "zmesh")

	(mesh_type, mesh_a, mesh_b) = guess_mesh_type(upf["radial_grid"],
												  upf["radial_grid_derivative"])
	if mesh_type == "unknown"
		raise(ExceptionError("Unknown mesh type"))
	end
	return upf["radial_grid_parameters"] = Dict("dx" => dx,
												"mesh" => mesh,
												"xmin" => xmin,
												"rmax" => rmax,
												"zmesh" => zmesh,
												"a" => mesh_a,
												"b" => mesh_b,
												"mesh_type" => mesh_type)
end

"""
	parse_nlcc_upf2!(io::IO, upf::Dict)

Parse non-linear core correction data from the `<PP_NLCC>` bock if present, storing them in
`upf["core_charge_density"]`.

``Z_c = \\int \\rho_c(r) r^2 dr d\\Omega``
"""
function parse_nlcc_upf2!(doc_root::EzXML.Node, upf::Dict)
	if upf["header"]["core_correction"]
		node = findfirst("PP_NLCC", doc_root)
		# Z_c = ∫(ρ_c(r) r^2 dr dΩ)
		upf["core_charge_density"] = parse.(Float64, split(strip(node.content)))
	end
end

"""
	parse_local_upf2!(io::IO, upf::Dict)

Parse the local potential from the `<PP_LOCAL>` block, storing it in
`upf["local_potential"]`. The local potential contains a long-range term
``-Z_\\text{valence} e^2 / r``.
"""
function parse_local_upf2!(doc_root::EzXML.Node, upf::Dict)
	node = findfirst("PP_LOCAL", doc_root)
	# Contains long range term -z_valence * e^2 / r
	return upf["local_potential"] = parse.(Float64, split(strip(node.content)))  # Ry
end

"""
	parse_beta_projectors_upf2!(io::IO, upf::Dict)

Parse the `<PP_BETA>` blocks in `<PP_NONLOCAL>`, storing them in a vector in
`upf["beta_projectors"]`. There are `upf["header"]["number_of_proj"]` blocks,
each with the following data:
- `label::String`: optional descriptive label
- `index::Int`: index of the projector, used for correlating with Dij
- `angular_momentum::Int`
- `cutoff_radius_index::Int`: number of elements read from file, all others are zero
- `radial_function::Vector{Float64}`: the beta projector, with length `cutoff_radius_index`
- `cutoff_radius::Float64`: always `0.`
- `ultrasoft_cutoff_radius::Float64`: always `0.`

If `upf["header"]["has_so"]`, spin-orbit data are parsed for each beta from
`PP_SPIN_ORBL/PP_RELBETA.[i]`, overwriting the data from `PP_NONLOCAL/PP_BETA.[i]`:
- `index::Int`: index of the projector, used for correlating with Dij
- `angular_momentum::Int`
- `total_angular_momentum::Int`

!!! Note
The units of the projectors are either ``\\text{Bohr}^{-1/2}`` or
``\\text{Ry}~\\text{Bohr}^{-1/2}``.
"""
function parse_beta_projectors_upf2!(doc_root::EzXML.Node, upf::Dict)
	beta_projectors = []
	for i in 1:upf["header"]["number_of_proj"]
		node = findfirst("PP_NONLOCAL/PP_BETA.$i", doc_root)
		beta = Dict()
		# Units are either Bohr^(-1/2) or Ry*Bohr^(-1/2)
		# The quantity is actually rᵢβ(rᵢ)
		beta["label"] = get_attr(String, node, "label")
		beta["angular_momentum"] = get_attr(Int, node, "angular_momentum")
		ir_cut = get_attr(Int, node, "cutoff_radius_index")
		beta["cutoff_radius_index"] = ir_cut
		beta["cutoff_radius"] = get_attr(Float64, node, "cutoff_radius")
		beta["index"] = get_attr(Int, node, "index")
		beta["radial_function"] = parse.(Float64, split(strip(node.content)))[1:ir_cut]
		if upf["header"]["has_so"]
			node_so = findfirst("PP_SPIN_ORB/PP_RELBETA.$i", doc_root)
			beta["index"] = get_attr(Int, node_so, "index")
			beta["angular_momentum"] = get_attr(Int, node_so, "lll")
			beta["total_angular_momentum"] = get_attr(Float64, node_so, "jjj")
		end
		push!(beta_projectors, beta)
	end
	return upf["beta_projectors"] = beta_projectors
end

"""
	parse_dij_upf2!(io::IO, upf::Dict)

Parse the `<PP_DIJ>` block, storing it in `upf["D_ion"]` as a symmetric matrix where
`D[i,j]` is the coupling coefficient between βᵢ and βⱼ (see `parse_beta_projectors_upf2!`
for how the indices of the beta projectors are stored).

!!! Note
The units of ``D_{ij}`` are either ``\\text{Ry}`` or ``\\text{Ry}^-1``, corresponding to the
units of the projectors.
"""
function parse_dij_upf2!(doc_root::EzXML.Node, upf::Dict)
	node = findfirst("PP_NONLOCAL/PP_DIJ", doc_root)
	Dij = parse.(Float64, split(strip(node.content)))
	Dij = reshape(Dij, upf["header"]["number_of_proj"], upf["header"]["number_of_proj"])
	return upf["D_ion"] = Dij  # either Ry or Ry^-1
end

# function parse_augmentation_upf2!(doc_root::EzXML.Node, upf::Dict)
# 	if !upf["header"]["is_ultrasoft"]
# 		augmentation = []
# 	else
# 		node = findfirst("PP_NONLOCAL/PP_AUGMENTATION", doc_root)
# 		q_with_l = get_attr(Bool, node, "q_with_l")
# 		if ismissing(q_with_l)
# 			throw(ErrorException("Parsing `q_with_l = T` is not supported."))
# 		end

# 		augmentation = []
# 		for i in 1:upf["header"]["number_of_proj"]
# 			li = upf["beta_projectors"][i]["angular_momentum"]
# 			for j in i:upf["header"]["number_of_proj"]
# 				lj = upf["beta_projectors"][j]["angular_momentum"]
# 				for l in abs(li - lj):(li + lj)
# 					if (li + lj + l) % 2 == 0
# 						node = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_QIJL.$i.$j.$l",
# 										 doc_root)
# 						Qij = Dict()
# 						Qij["radial_function"] = parse.(Float64, split(strip(node.content)))
# 						Qij["i"] = i
# 						Qij["j"] = j
# 						Qij["angular_momentum"] = get_attr(Int, node, "angular_momentum")
# 						push!(augmentation, Qij)
# 					end
# 				end
# 			end
# 		end
# 	end
# 	return upf["augmentation"] = augmentation
# end

# function parse_paw_upf2!(doc_root::EzXML.Node, upf::Dict)
# 	if !(lowercase(upf["header"]["pseudo_type"]) == "paw")
# 		paw_data = Dict()
# 	else
# 		node = findfirst("PP_NONLOCAL/PP_AUGMENTATION", doc_root)
# 		upf["header"]["cutoff_radius_index"] = get_attr(Int, node, "cutoff_r_index")

# 		paw_data = Dict()

# 		node_q = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_Q", doc_root)
# 		paw_data["aug_integrals"] = parse.(Float64, split(strip(node_q.content)))

# 		node_mp = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_MULTIPOLES", doc_root)
# 		paw_data["aug_multipoles"] = parse.(Float64, split(strip(node_mp.content)))

# 		paw_data["ae_wfc"] = []
# 		for i in 1:upf["header"]["number_of_proj"]
# 			wfc = Dict()
# 			node_wfc = findfirst("PP_FULL_WFC/PP_AEWFC.$i", doc_root)
# 			wfc["radial_function"] = parse.(Float64, split(strip(node_wfc.content)))
# 			wfc["angular_momentum"] = get_attr(Int, node_wfc, "l")
# 			wfc["label"] = get_attr(String, node_wfc, "label")
# 			wfc["index"] = get_attr(Int, node_wfc, "index")
# 			push!(paw_data["ae_wfc"], wfc)
# 		end

# 		paw_data["ps_wfc"] = []
# 		for i in 1:upf["header"]["number_of_proj"]
# 			wfc = Dict()
# 			node_wfc = findfirst("PP_FULL_WFC/PP_PSWFC.$i", doc_root)
# 			wfc["radial_function"] = parse.(Float64, split(strip(node_wfc.content)))
# 			wfc["angular_momentum"] = get_attr(Int, node_wfc, "l")
# 			# wfc["label"] = get_attr(String, node_wfc, "label")
# 			# wfc["index"] = get_attr(Int, node_wfc, "index")
# 			push!(paw_data["ps_wfc"], wfc)
# 		end

# 		node_paw = findfirst("PP_PAW", doc_root)
# 		paw_core_energy = get_attr(Float64, node_paw, "core_energy"; default=missing)
# 		if ismissing(paw_core_energy)
# 			# @warn "`PP_PAW` has no `core_energy` set"
# 		else
# 			upf["header"]["paw_core_energy"] = paw_core_energy  # Ry
# 		end

# 		node_occ = findfirst("PP_PAW/PP_OCCUPATIONS", doc_root)
# 		paw_data["occupations"] = parse.(Float64, split(strip(node_occ.content)))

# 		node_ae_nlcc = findfirst("PP_PAW/PP_AE_NLCC", doc_root)
# 		paw_data["ae_core_charge_density"] = parse.(Float64,
# 													split(strip(node_ae_nlcc.content)))

# 		node_ae_vloc = findfirst("PP_PAW/PP_AE_VLOC", doc_root)
# 		paw_data["ae_local_potential"] = parse.(Float64, split(strip(node_ae_vloc.content)))  # Ry
# 	end
# 	return upf["paw_data"] = paw_data
# end

"""
	parse_pswfc_upf2!(io::IO, upf::Dict)

Parse the pseudo-atomic wavefunctions in the `<PP_PSWFC>` block, storing them in a vector in
`upf["atomic_wave_functions"]`. There are `upf["header"]["number_of_wfc"]` blocks,
each with the following data:
- `label::String`: optional descriptive label, e.g. "2S"
- `angular_momentum::Int`
- `occupation::Float`
- `index::Int`
- `pseudo_energy::Float`
- `radial_function::Vector{Float64}`: the pseudo-atomic wavefunction on the full radial mesh

If `upf["header"]["has_so"]`, spin-orbit data are parsed for each pseudo-atomic orbital from
`PP_SPIN_ORBL/PP_RELWFC.[i]`, overwriting the data from `PP_PSWFC/PP_CHI.[i]`:
- `index::Int`
- `angular_momentum::Int`
- `total_angular_momentum::Int`
- `principal_quantum_number::Int`
"""
function parse_pswfc_upf2!(doc_root::EzXML.Node, upf::Dict)
	atomic_wave_functions = []
	for i in 1:upf["header"]["number_of_wfc"]
		node = findfirst("PP_PSWFC/PP_CHI.$i", doc_root)
		wfc = Dict()
		wfc["label"] = get_attr(String, node, "label")
		wfc["angular_momentum"] = get_attr(Int, node, "l")
		wfc["occupation"] = get_attr(Float64, node, "occupation")
		wfc["pseudo_energy"] = get_attr(Float64, node, "pseudo_energy")  # Ry
		wfc["index"] = get_attr(Int, node, "index")
		wfc["radial_function"] = parse.(Float64, split(strip(node.content)))
		if upf["header"]["has_so"]
			node_so = findfirst("PP_SPIN_ORB/PP_RELWFC.$i", doc_root)
			wfc["index"] = get_attr(Int, node_so, "index")
			wfc["angular_momentum"] = get_attr(Int, node_so, "lchi")
			wfc["total_angular_momentum"] = get_attr(Float64, node_so, "jchi")
			wfc["principal_quantum_number"] = get_attr(Int, node_so, "nn")
		end
		push!(atomic_wave_functions, wfc)
	end
	return upf["atomic_wave_functions"] = atomic_wave_functions
end

"""
	parse_rhoatom_upf2!(io::IO, upf::Dict)

Parse the total pseudo-atomic charge density from the `<PP_RHOATOM>` block, storing the data
in `upf["total_charge_density"]`.

!!! Note
There is _no_ ``4π`` prefactor!
"""
function parse_rhoatom_upf2!(doc_root::EzXML.Node, upf::Dict)
	node = findfirst("PP_RHOATOM", doc_root)
	return upf["total_charge_density"] = parse.(Float64, split(strip(node.content)))
end

# function parse_spin_orbit_upf2!(doc_root::EzXML.Node, upf::Dict)
# 	if upf["header"]["has_so"]
# 		for i in 1:upf["header"]["number_of_proj"]
# 			node = findfirst("PP_SPIN_ORB/PP_RELBETA.$i", doc_root)
# 			upf["beta_projectors"][i]["angular_momentum"] = get_attr(Float64, node, "lll")
# 			upf["beta_projectors"][i]["total_angular_momentum"] = get_attr(Float64, node,
# 																		   "jjj")
# 		end

# 		for i in 1:upf["header"]["number_of_wfc"]
# 			node = findfirst("PP_SPIN_ORB/PP_RELWFC.$i", doc_root)
# 			upf["paw_data"]["ae_wfc"]["ae_wfc_rel"] = parse(Float64, strip(node.content))
# 			upf["paw_data"]["ae_wfc"]["total_angular_momentum"] = get_attr(Float64, node,
# 																		   "jchi")
# 			upf["paw_data"]["ps_wfc"]["total_angular_momentum"] = get_attr(Float64, node,
# 																		   "jchi")
# 		end
# 	end
# end

"""
	parse_upf2(doc_root::EzXML.Node)

Parse a UPF v2 (with schema) file.

All quantities are in Rydberg units:
- e² = 2
- m = 1 / 2
- ħ = 1
- Lengths in Bohr (0.529177 Å)
- Energies in Ry (13.6058 eV)
- Potentials multiplied by e to give units of energy

!!! Note
PAW and ultrasoft potentials are not supported because parsing of `<PP_AUGMENTATION>` and
`<PP_PAW>` are not fully implemented.
"""
function parse_upf2(doc_root::EzXML.Node)
	upf = Dict()

	parse_header_upf2!(doc_root, upf)
	if upf["header"]["pseudo_type"] == "PAW"
		@warn "PAW in UPF v2 is not implemented."
	elseif upf["header"]["pseudo_type"] == "US"
		@warn "Ultrasoft in UPF v2 is not implemented."
	end
	parse_radial_grid_upf2!(doc_root, upf)
	parse_nlcc_upf2!(doc_root, upf)
	parse_local_upf2!(doc_root, upf)
	parse_beta_projectors_upf2!(doc_root, upf)
	parse_dij_upf2!(doc_root, upf)
	# parse_augmentation_upf2!(doc_root, upf)
	# parse_paw_upf2!(doc_root, upf)
	parse_pswfc_upf2!(doc_root, upf)
	parse_rhoatom_upf2!(doc_root, upf)
	# parse_spin_orbit_upf2!(doc_root, upf)

	return upf
end
