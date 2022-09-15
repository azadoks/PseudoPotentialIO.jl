using LinearAlgebra

function parse_header(io)
    header = Dict()
    
    # line 1: title
    header["title"] = readline(io)              # title (unused)
    
    # line 2: atomic number, pseudoion charge, date
    s = split(readline(io))
    header["element"] = parse(Float64, s[1])    # zatom
    header["z_valence"] = parse(Float64, s[2])  # zion
    header["date"] = parse(Int, s[3])           # pspd (unused)
    
    # line 3
    s = split(readline(io))
    header["format_version"] = parse(Int, s[1]) # pspcod == 8
    header["xc"] = parse(Int, s[2])             # pspxc
    header["l_max"] = parse(Int, s[3])          # lmax
    header["l_local"] = parse(Int, s[4])        # lloc
    header["mesh_size"] = parse(Int, s[4])      # mmax
    header["r2well"] = parse(Int, s[4])         # r2well (unused)
    @assert header["format_version"] == 8

    # line 4
    s = split(readline(io))
    header["rchrg"] = parse(Float64, s[1])       # rchrg
    header["fchrg"] = parse(Float64, s[2])       # fchrg
    header["qchrg"] = parse(Float64, s[3])       # qchrg (unused)
    header["core_correction"] = header["fchrg"] > 0.0

    # line 5: number of scalar-relativistic non-local
    # projectors for each angular momentum (l = 0:l_max)
    s = split(readline(io))
    header["number_of_proj"] = [parse(Int, s[i+1]) for i = 0:header["l_max"]]
    if header["l_local"] <= header["l_max"]
        @assert header["number_of_proj"][header["l_local"]+1] == 0
    end

    # line 6: data extension information
    s = split(readline(io))
    header["extension_switch"] = []
    for s_i in s
        try
            push!(header["extension_switch"], parse(Int, s_i))
        catch
        end
    end
    header["has_so"] = header["extension_switch"][1] == 2 ? true : false

    if header["has_so"]
        # line 7: number of projectors for each spin-orbit
        # non-local projectors for each angular momentum (l = 1:l_max)
        header["number_of_proj_so"] = [parse(Int, s[i+1]) for i = 1:header["l_max"]]
    else
        header["number_of_proj_so"] = []
    end

    return header
end

function parse_beta_projector(io, psp)
    header_line = split(readline(io))
    l = parse(Int, header_line[1])
    n_proj_l = header["number_of_proj"][l + 1]
    ekb = [parse(Float64, header_line[i]) for i = 2:(2 + n_proj_l)]
    
    radial_grid = Vector{Float64}(undef, psp["header"]["mesh_size"])
    betas = [Vector{Float64}(undef, psp["header"]["mesh_size"]) for i = 1:(n_proj_l)]
    for i = 1:psp["header"]["mesh_size"]
        s = split(readline(io))
        radial_grid[i] = parse(Float64, s[2])
        for j = 1:n_proj_l
            betas[j][i] = parse(Float64, s[2+j])
        end
    end

    return Dict(
        "angular_momentum" => l,
        "radial_grid" => radial_grid,
        "radial_functions" => betas,
        "ekb" => ekb
    )
end

function parse_local(io, psp)
    header_line = split(readline(io))
    l = parse(Int, header_line[1])

    radial_grid = Vector{Float64}(undef, psp["header"]["mesh_size"])
    v_local = Vector{Float64}(undef, psp["header"]["mesh_size"])
    for i = 1:psp["header"]["mesh_size"]
        s = split(readline(io))
        radial_grid[i] = parse(Float64, s[2])
        v_local[i] = parse(Float64, s[3])
    end
    return Dict(
        "angular_momentum" => l,
        "radial_grid" => radial_grid,
        "local_potential" => v_local
    )
end

function parse_betas_local(io, psp)
    beta_blocks = []
    v_local_block = Dict()
    for i = 1:(psp["header"]["l_max"] + 1)
        # Record the position at the start of the block so we can
        # read in the first line and go back
        block_head = position(io)
        # Read the block header
        block_header_line = split(readline(io))
        # Go back to the start of the block
        seek(io, block_head)
        # Parse the block's `l`
        block_l = parse(Int, block_header_line[1])

        if block_l == psp["header"]["l_local"]
            v_local_block = parse_local(io, psp)
        else
            beta_block = parse_beta_projector(io, psp)
            push!(beta_blocks, beta_block)
        end
    end

    beta_projectors = [
        Dict(
            "radial_function" => beta["radial_function"],
            "angular_momentum" => beta["angular_momentum"],
        ) for beta in beta_blocks
    ]

    Dij = []
    for beta_block in beta_blocks
        for dii in beta_block["ekb"]
            push!(Dij, dii)
        end
    end
    Dij = collect(Diagonal(Dij))

    return Dict(
        "radial_grid" => v_local_block["radial_grid"],
        "local_potential" => v_local_block["local_potential"],
        "beta_projectors" => beta_projectors,
        "D_ion" => Dij
    )
end

function parse_nlcc(io, psp)
    mesh_size = psp["header"]["mesh_size"]
    radial_grid = Vector{Float64}(undef, mesh_size)
    rho = Vector{Float64}(undef, mesh_size)
    drho = Vector{Float64}(undef, mesh_size)
    d2rho = Vector{Float64}(undef, mesh_size)
    d3rho = Vector{Float64}(undef, mesh_size)
    d4rho = Vector{Float64}(undef, mesh_size)
    for i = 1:mesh_size
        s = split(readline(io))
        radial_grid[i] = parse(Float64, s[2])
        rho[i] = parse(Float64, s[2])
        drho[i] = parse(Float64, s[3])
        d2rho[i] = parse(Float64, s[4])
        d3rho[i] = parse(Float64, s[5])
        d4rho[i] = parse(Float64, s[6])
    end
    return Dict(
        "radial_grid" => radial_grid,
        "rho" => rho,
        "drho" => drho,
        "d2rho" => d2rho,
        "d3rho" => d3rho,
        "d4rho" => d4rho
    )
end

function parse_so(io, psp)
end

function parse_psp8(io)
    psp = Dict()

    psp["header"] = parse_header(io)
    for (key, value) in parse_betas_local(io, psp)
        psp[key] = value
    end

    if psp["header"]["has_so"]
        so_blocks = []
        for i = 1:(psp["header"]["l_max"])
            so_block = parse_so(io, psp)
            push!(so_blocks, so_block)
        end
    else
        psp["atomic_wave_functions"] = []
    end

    if psp["header"]["core_correction"]
        nlcc = parse_nlcc(io, psp)
        psp["core_charge_density"] = nlcc["core_charge_density"]
    else
        psp["core_charge_density"] = []
    end
end