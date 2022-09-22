module UPF1

export parse_upf1

function read_until(io, tag::AbstractString)
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

function read_mesh_data(T::Type, io::IO, n::Integer)
    mesh_data = T[]
    while length(mesh_data) != n
        line_data = parse.(T, split(readline(io)))
        append!(mesh_data, line_data)
    end
    return mesh_data
end;

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

function parse_header!(io::IO, upf::Dict)
    header = Dict()

    read_until(io, "<PP_HEADER>")

    header["format_version"] = 1

    s = split(readline(io))
    header["version"] = parse(Int, s[1])        # Version Number

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

function parse_radial_grid!(io::IO, upf::Dict)
    read_until(io, "<PP_R>")
    upf["radial_grid"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])

    read_until(io, "<PP_RAB>")
    upf["radial_grid_derivative"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
end

function parse_nlcc!(io::IO, upf::Dict)
    if !upf["header"]["core_correction"]
        rho = missing
    else
        read_until(io, "<PP_NLCC>")
        rho = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
    end
    upf["core_charge_density"] = rho
end

function parse_local!(io::IO, upf::Dict)
    read_until(io, "<PP_LOCAL>")
    upf["local_potential"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"]) ./ 2  # Ry -> Ha
end

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
        beta["ultrasoft_cutoff_radius"] = 0.
        
        push!(beta_projectors, beta)
    end
    
    upf["beta_projectors"] = beta_projectors
end

function parse_dij!(io::IO, upf::Dict)
    read_until(io, "<PP_DIJ>")
    Dij = zeros(Float64, upf["header"]["number_of_proj"], upf["header"]["number_of_proj"])
    s = split(readline(io))
    nd = parse(Int, s[1])  # Number of non-zero Dij components
    for k = 1:nd
        s = split(readline(io))
        i, j = parse.(Int, s[1:2])
        Dij[i,j] = parse(Float64, s[3]) / 2  # Ry -> Ha
        Dij[j,i] = Dij[i,j]
    end
    
    upf["D_ion"] = Dij
end

function parse_augmentation!(io::IO, upf::Dict)
    if !(uppercase(upf["header"]["pseudo_type"]) != "US" || uppercase(upf["header"]["pseudo_type"] == "PAW"))
        upf["augmentation"] = missing
        return
    end
    
    read_until(io, "<PP_QIJ>")
    # If num_q_coef is non-zero, Qij inside R_inner (parsed below,
    # one value for each augmentation) are computed using the q_coeffs
    num_q_coeff = parse(Int, split(readline(io))[1])
    if num_q_coeff == 0
        upf["augmentation"] = missing
        return
    end

    augmentation = []

    R_inner = Float64[]  # Inner radial cutoff values for each augmentation
    read_until(io, "<PP_RINNER>")
    for i = 1:(2 * upf["header"]["l_max"] + 1)
        s = split(readline(io))
        push!(R_inner, parse(Float64, s[2]))
    end
    read_until(io, "</PP_RINNER>")

    for i = 1:upf["header"]["number_of_proj"]
        li = upf["beta_projectors"][i]["angular_momentum"]
        for j = i:upf["header"]["number_of_proj"]
            lj = upf["beta_projectors"][j]["angular_momentum"]

            s = split(readline(io))
            file_i = parse(Int, s[1])
            file_j = parse(Int, s[2])
            file_lj = parse(Int, s[3])

            s = split(readline(io))
            q_int = parse(Float64, s[1])

            qij = read_mesh_data(Float64, io, upf["header"]["mesh_size"])

            read_until(io, "<PP_QFCOEF>")
            q_coeffs = read_mesh_data(Float64, io, num_q_coeff * (2 * upf["header"]["l_max"] + 1))
            q_coeffs = reshape(q_coeffs, num_q_coeff, 2 * upf["header"]["l_max"] + 1)
            read_until(io, "</PP_QFCOEF>")

            for l = abs(li - lj):(li + lj)
                if (li + lj + l) % 2 == 0
                    qij_fixed = copy(qij)

                    for ir = 1:upf["header"]["mesh_size"]
                        x = upf["radial_grid"][ir]
                        if x < R_inner[l+1]
                            qij_fixed[ir] = sum(n -> q_coeffs[n,l+1] * x^(2 * (n-1)), 1:num_q_coeff)
                            qij_fixed[ir] *= x^(l + 2)
                        end
                    end
                    
                    push!(augmentation, Dict(
                        "radial_function" => qij_fixed,
                        "i" => i,
                        "j" => j,
                        "angular_momentum" => l,
                    ))
                end
            end
        end
    end
    upf["augmentation"] = augmentation
end

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

function parse_rhoatom!(io::IO, upf::Dict)
    read_until(io, "<PP_RHOATOM>")
    upf["total_charge_density"] = read_mesh_data(Float64, io, upf["header"]["mesh_size"])
end

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

function parse_upf1(io::IO)
    upf = Dict()

    parse_header!(io, upf)
    if upf["header"]["is_paw"]
        @warn "PAW is not implemented."
    elseif upf["header"]["is_ultrasoft"]
        @warn "Ultrasoft is not implemented."
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