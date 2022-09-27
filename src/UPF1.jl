module UPF1

export parse_upf1

"""
    read_until(io, tag::AbstractString)

Read until a line contains `tag`.
"""
function read_until(io::IO, tag::AbstractString)
    while true
        line = readline(io, keep=true)
        if isempty(line)
            throw(EOFError())
        end
        if occursin(tag, line)
            return
        end
    end
end;

"""
    read_mesh_data(T::Type, io::IO, n::Integer)

Read whitespace-separated data from `io` as type `T` until `n` elements are read.
"""
function read_mesh_data(T::Type, io::IO, n::Integer)
    mesh_data = T[]
    while length(mesh_data) != n
        line_data = parse.(T, split(readline(io)))
        append!(mesh_data, line_data)
    end
    return mesh_data
end;

"""
    has_tag(io::IO, tag::AbstractString)

Check if `tag` occurs in `io` after the current position in `io`.
"""
function has_tag(io::IO, tag::AbstractString)::Bool
    pos = position(io)
    has_tag = false
    try
        read_until(io, tag)
        has_tag = true
        seek(io, pos)
    catch
        has_tag = false
        seek(io, pos)
    end
    return has_tag
end

"""
    parse_header!(io::IO, upf::Dict)

Parse header (`PP_HEADER`) data, storing it in `upf["header"]::Dict` with the following contents:
- `format_version::Int`: always 1
- `version::Int`: version number of pseudopotential (not format!)
- `element::String`: elemental symbol
- `pseudo_type::String`: `NC` (norm-conserving), `US` (ultrasoft), or `PAW` (plane-augmented wave)
- `is_ultrasoft::Bool`
- `is_paw::Bool`
- `is_coulomb::Bool`: fake Coulomb potential for all-electron calculations
- `has_so::Bool`: has spin-orbit coupling
- `has_wfc::Bool`: always false
- `has_gipaw::Bool`: has GIPAW reconstruction data
- `core_correction::Bool`: non-linear core-correction
- `z_valence::Float64`: pseudo-ion charge
- `total_psenergy::Float64`: total energy of the pseudo-ion
- `ecutwfc::Float64`: recommended plane-wave energy cutoff
- `ecutrho::Float64`: recommended charge-density energy cutoff
- `l_max::Int`: maximum angular momentum channel of the Kleinman-Bylander projectors
- `mesh_size::Int`: number of points in the radial mesh
- `number_of_wfc::Int`: number of pseudo-atomic wavefunctions
- `number_of_proj::Int`: number of Kleinman-Bylander projectors

- `q_with_l::Bool`: Qij is dependent on angular momentum

!!! Note
The "Wavefunctions" section of the header is not parsed; instead the label, angular momentum,
and occupations are parsed with the wavefunctions in the "PP_PSWFC" block.
"""
function parse_header!(io::IO, upf::Dict)
    header = Dict()

    read_until(io, "<PP_HEADER>")

    header["format_version"] = 1

    s = split(readline(io))
    header["version"] = parse(Int, s[1])        # Version number

    s = split(readline(io))
    header["element"] = strip(s[1])             # Element

    s = split(readline(io))
    header["pseudo_type"] = s[1]                # Type of pseudopotential (NC|US|PAW)
    header["is_ultrasoft"] = s[1] == "US"
    header["is_paw"] = s[1] == "PAW"
    header["is_coulomb"] = s[1] == "1/r"

    s = split(readline(io))
    header["core_correction"] = occursin("T", uppercase(s[1])) ? true : false

    s = split(readline(io))
    # header["dft"] = s[begin:20]                 # Exchange-correlation functional

    s = split(readline(io))
    header["z_valence"] = parse(Float64, s[1])  # Z valence

    s = split(readline(io))
    header["total_psenergy"] = parse(Float64, s[1])     # Total energy

    s = split(readline(io))
    header["ecutwfc"] = parse(Float64, s[1])    # Suggested cutoffs
    header["ecutrho"] = parse(Float64, s[2])

    s = split(readline(io))
    header["l_max"] = parse(Int, s[1])          # Maximum angular momentum

    s = split(readline(io))
    header["mesh_size"] = parse(Int, s[1])      # Number of points in radial mesh

    s = split(readline(io))
    header["number_of_wfc"] = parse(Int, s[1])  # Number of wavefunctions
    header["number_of_proj"] = parse(Int, s[2]) # Number of projectors

    #! We don't parse the "Wavefunctions" description block

    header["has_so"] = has_tag(io, "<PP_ADDINFO>")
    header["has_gipaw"] = has_tag(io, "<PP_GIPAW_RECONSTRUCTION_DATA>")
    header["q_with_l"] = has_tag(io, "<PP_QIJ_WITH_L>")
    header["has_wfc"] = false

    upf["header"] = header
end

"""
    parse_radial_grid!(io::IO, upf::Dict)

Parse radial grid data (`<PP_R>`) and integration factors (`<PP_RAB>`) from the `<PP_MESH>` block,
storing them in `upf["radial_grid"]` and `upf["radial_grid_derivative"]` respectively.

The radial grid is one of the following:
``e^{x_\\text{min}} e^{(i - 1)dx} / Z_\\text{mesh}``
``e^{x_\\text{min}} (e^{(i - 1)dx} - 1) / Z_\\text{mesh}``

The radial grid derivative is the factor for discrete integration:
``\\int f(r) dr = \\sum_{i=1}^{N} f(i) r_{ab}(i)``  
"""
function parse_radial_grid!(io::IO, upf::Dict)
    read_until(io, "<PP_R>")
    upf["radial_grid"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])

    read_until(io, "<PP_RAB>")
    upf["radial_grid_derivative"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
end

"""
    parse_nlcc!(io::IO, upf::Dict)

Parse non-linear core correction data from the `<PP_NLCC>` bock if present, storing them in
`upf["core_charge_density"]`.

``Z_c = \\int \\rho_c(r) r^2 dr d\\Omega``
"""
function parse_nlcc!(io::IO, upf::Dict)
    if upf["header"]["core_correction"]
        read_until(io, "<PP_NLCC>")
        upf["core_charge_density"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
    end
end

"""
    parse_local!(io::IO, upf::Dict)

Parse the local potential from the `<PP_LOCAL>` block, storing it in `upf["local_potential"]`.
"""
function parse_local!(io::IO, upf::Dict)
    read_until(io, "<PP_LOCAL>")
    upf["local_potential"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
end

"""
    parse_beta_projectors!(io::IO, upf::Dict)

Parse the `<PP_BETA>` blocks in `<PP_NONLOCAL>`, storing them in a vector in
`upf["beta_projectors"]`. There are `upf["header"]["number_of_proj"]` blocks,
each with the following data:
- `label::String`: optional descriptive label
- `index::Int`: index of the projector, used for correlating with Dij
- `angular_momentum::Int`
- `cutoff_radius_index::Int`: number of elements read from file, all others are zero
- `radial_function::Vector{Float64}`: the beta projector, with length `cutoff_radius_index`
- `cutoff_radius::Float64`: always `0.`

!!! Note
The units of the projectors are either ``\\text{Bohr}^{-1/2}`` or
``\\text{Ry}~\\text{Bohr}^{-1/2}``.
"""
function parse_beta_projectors!(io::IO, upf::Dict)
    read_until(io, "<PP_NONLOCAL>")
    beta_projectors = []
    for i = 1:upf["header"]["number_of_proj"]
        read_until(io, "<PP_BETA>")
        
        beta = Dict()
        beta["label"] = ""
        
        s = split(readline(io))
        beta["index"] = parse(Int, s[1])
        beta["angular_momentum"] = parse(Int, s[2])

        s = readline(io)
        beta["cutoff_radius_index"] = parse(Int, s)

        beta["radial_function"] = read_mesh_data(Float64, io, beta["cutoff_radius_index"])
        beta["cutoff_radius"] = 0.
        
        push!(beta_projectors, beta)
    end
    
    upf["beta_projectors"] = beta_projectors
end

"""
    parse_dij!(io::IO, upf::Dict)

Parse the `<PP_DIJ>` block, storing it in `upf["D_ion"]` as a symmetric matrix where `D[i,j]` is
the coupling coefficient between βᵢ and βⱼ (see `parse_beta_projectors!` for how the indices of
the beta projectors are stored).

!!! Note
The units of ``D_{ij}`` are either ``\\text{Ry}`` or ``\\text{Ry}^-1``, corresponding to
the units of the projectors.
"""
function parse_dij!(io::IO, upf::Dict)
    read_until(io, "<PP_DIJ>")
    Dij = zeros(Float64, upf["header"]["number_of_proj"], upf["header"]["number_of_proj"])
    s = split(readline(io))
    nd = parse(Int, s[1])  # Number of non-zero Dij components
    for _ = 1:nd
        s = split(readline(io))
        i, j = parse.(Int, s[1:2])
        Dij[i,j] = parse(Float64, s[3])
        Dij[j,i] = Dij[i,j]
    end
    
    upf["D_ion"] = Dij
end

# function parse_augmentation!(io::IO, upf::Dict)
#     if !(
#         uppercase(upf["header"]["pseudo_type"]) != "US" ||
#         uppercase(upf["header"]["pseudo_type"] == "PAW")
#     )
#         upf["augmentation"] = missing
#         return
#     end
    
#     read_until(io, "<PP_QIJ>")
#     # If num_q_coef is non-zero, Qij inside R_inner (parsed below,
#     # one value for each augmentation) are computed using the q_coeffs
#     num_q_coeff = parse(Int, split(readline(io))[1])
#     if num_q_coeff == 0
#         upf["augmentation"] = missing
#         return
#     end

#     augmentation = []

#     R_inner = Float64[]  # Inner radial cutoff values for each augmentation
#     read_until(io, "<PP_RINNER>")
#     for i = 1:(2 * upf["header"]["l_max"] + 1)
#         s = split(readline(io))
#         push!(R_inner, parse(Float64, s[2]))
#     end
#     read_until(io, "</PP_RINNER>")

#     for i = 1:upf["header"]["number_of_proj"]
#         li = upf["beta_projectors"][i]["angular_momentum"]
#         for j = i:upf["header"]["number_of_proj"]
#             lj = upf["beta_projectors"][j]["angular_momentum"]

#             s = split(readline(io))
#             file_i = parse(Int, s[1])
#             file_j = parse(Int, s[2])
#             file_lj = parse(Int, s[3])

#             s = split(readline(io))
#             q_int = parse(Float64, s[1])

#             qij = read_mesh_data(Float64, io, upf["header"]["mesh_size"])

#             read_until(io, "<PP_QFCOEF>")
#             q_coeffs = read_mesh_data(Float64, io, num_q_coeff * (2 * upf["header"]["l_max"] + 1))
#             q_coeffs = reshape(q_coeffs, num_q_coeff, 2 * upf["header"]["l_max"] + 1)
#             read_until(io, "</PP_QFCOEF>")

#             for l = abs(li - lj):(li + lj)
#                 if (li + lj + l) % 2 == 0
#                     qij_fixed = copy(qij)

#                     for ir = 1:upf["header"]["mesh_size"]
#                         x = upf["radial_grid"][ir]
#                         if x < R_inner[l+1]
#                             qij_fixed[ir] = sum(n -> q_coeffs[n,l+1] * x^(2 * (n-1)), 1:num_q_coeff)
#                             qij_fixed[ir] *= x^(l + 2)
#                         end
#                     end
                    
#                     push!(augmentation, Dict(
#                         "radial_function" => qij_fixed,
#                         "i" => i,
#                         "j" => j,
#                         "angular_momentum" => l,
#                     ))
#                 end
#             end
#         end
#     end
#     upf["augmentation"] = augmentation
# end

"""
    parse_pswfc!(io::IO, upf::Dict)

Parse the pseudo-atomic wavefunctions in the `<PP_PSWFC>` block, storing them in a vector in
`upf["atomic_wave_functions"]`. There are `upf["header"]["number_of_wfc"]` blocks,
each with the following data:
- `label::String`: optional descriptive label, e.g. "2S"
- `angular_momentum::Int`
- `occupation::Float`
- `radial_function::Vector{Float64}`: the pseudo-atomic wavefunction on the full radial mesh
"""
function parse_pswfc!(io::IO, upf::Dict)
    atomic_wave_functions = []
    read_until(io, "<PP_PSWFC>")
    for i = 1:upf["header"]["number_of_wfc"]
        wfc = Dict()
        s = split(readline(io))
        wfc["label"] = s[1]
        wfc["angular_momentum"] = parse(Int, s[2])
        wfc["occupation"] = parse(Float64, s[3])
        wfc["radial_function"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
        push!(atomic_wave_functions, wfc)
    end
    upf["atomic_wave_functions"] = atomic_wave_functions
end

"""
    parse_rhoatom!(io::IO, upf::Dict)

Parse the total pseudo-atomic charge density from the `<PP_RHOATOM>` block, storing the data in
`upf["total_charge_density"]`.
"""
function parse_rhoatom!(io::IO, upf::Dict)
    read_until(io, "<PP_RHOATOM>")
    upf["total_charge_density"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
end

"""
    parse_addinfo!(io::IO, upf::Dict)

If `upf["header"]["has_so"]`, parse spin-orbit coupling data, which are stored in the
`<PP_ADDINFO>` block. For each pseudo-atomic wavefunction, the `label`, `angular_momentum`, and
`occupation` are overwritten, and new keys `principal_quantum_number` and `total_angular_momentum`
are stored. For each Kleinman-Bylander projector, `angular_momentum` is overwritten, and a new
`total_angular_momentum` is stored.
"""
function parse_addinfo!(io::IO, upf::Dict)
    if upf["header"]["has_so"]
        read_until(io, "<PP_ADDINFO>")
        for i = 1:upf["header"]["number_of_wfc"]
            s = split(readline(io))
            label = strip(s[1]) 
            n = parse(Int, s[2])
            l = parse(Int, s[3])
            j = parse(Float64, s[4]) 
            occ = parse(Float64, s[5])
            @assert (abs(j - l - 0.5) ≈ 0) || (abs(j - l + 0.5) ≈ 0)
            upf["atomic_wave_functions"][i]["label"] = label
            upf["atomic_wave_functions"][i]["principal_quantum_number"] = n
            upf["atomic_wave_functions"][i]["total_angular_momentum"] = j
            upf["atomic_wave_functions"][i]["angular_momentum"] = l
            upf["atomic_wave_functions"][i]["occupation"] = occ
        end

        for i = 1:upf["header"]["number_of_proj"]
            s = split(readline(io))
            l = parse(Int, s[1])
            j = parse(Float64, s[2])
            @assert ((abs(j - l - 0.5) ≈ 0) || (abs(j - l + 0.5) ≈ 0))
            upf["beta_projectors"][i]["angular_momentum"] = l
            upf["beta_projectors"][i]["total_angular_momentum"] = j
        end

        s = split(readline(io))  # xmin, rmax, zmesh, dx
    end
end

"""
    parse_upf1(io::IO)

Parse an old UPF (v1) file.

All quantities are in Rydberg units:
- e² = 2
- m = 1 / 2
- ħ = 1
- Lengths in Bohr (0.529177 Å)
- Energies in Ry (13.6058 eV)
- Potentials multiplied by e to give units of energy

!!! Note
PAW and ultrasoft potentials are not supported because parsing of `<PP_AUGMENTATION>` is not
fully implemented.
"""
function parse_upf1(io::IO)
    upf = Dict()

    parse_header!(io, upf)
    if upf["header"]["is_paw"]
        @warn "PAW in UPF v1 is not implemented."
    elseif upf["header"]["is_ultrasoft"]
        @warn "Ultrasoft in UPF v1 is not implemented."
    end
    parse_radial_grid!(io, upf)
    parse_nlcc!(io, upf)
    parse_local!(io, upf)
    parse_beta_projectors!(io, upf)
    parse_dij!(io, upf)
    # parse_augmentation!(io, upf)
    parse_pswfc!(io, upf)
    parse_rhoatom!(io, upf)
    parse_addinfo!(io, upf)
    return upf
end
end