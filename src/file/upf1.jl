function read_until(io::IO, tag::AbstractString)
	while true
		line = readline(io; keep=true)
		if isempty(line)
			throw(EOFError())
		end
		if occursin(tag, line)
			return
		end
	end
end

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

function read_mesh_data(T::Type, io::IO, n::Integer)
    mesh_data = T[]
    while length(mesh_data) != n
        line_data = parse.(T, split(readline(io)))
        append!(mesh_data, line_data)
    end
    return mesh_data
end

function upf1_parse_info(io::IO)
    pos = position(io)
    info = String[]
    read_until(io, "<PP_INFO>")
    line = string(strip(readline(io)))
    push!(info, line)
    while !occursin("</PP_INFO>", line)
        line = string(strip(readline(io)))
        push!(info, line)
    end
    seek(io, pos)
    return join(info, '\n')
end

function upf1_parse_header(io::IO)
    read_until(io, "<PP_HEADER>")

    s = split(readline(io))
    _ = parse(Int, s[1])  # version

    s = split(readline(io))
    element = strip(s[1])

    s = split(readline(io))
    pseudo_type = s[1]
    is_ultrasoft = pseudo_type in ("PAW", "US", "USPP")
    is_paw = pseudo_type == "PAW"
    is_coulomb = pseudo_type == "1/r"
    pseudo_type == "SL" && error("Semilocal pseudopotentials are not supported")

    s = split(readline(io))
    core_correction = occursin("T", uppercase(s[1])) ? true : false

    s = strip.(split(strip(readline(io))))
    functional = join(filter(s_i -> length(s_i) <= 6 , s), ' ')

    s = split(readline(io))
    z_valence = parse(Float64, s[1])

    s = split(readline(io))
    total_psenergy = parse(Float64, s[1])

    s = split(readline(io))
    wfc_cutoff = parse(Float64, s[1])
    rho_cutoff = parse(Float64, s[2])

    s = split(readline(io))
    l_max = parse(Int, s[1])

    s = split(readline(io))
    mesh_size = parse(Int, s[1])

    s = split(readline(io))
    number_of_wfc = parse(Int, s[1])
    number_of_proj = parse(Int, s[2])

    generated = nothing
    author = nothing
    date = nothing
    comment = nothing
    paw_as_gipaw = nothing
    l_max_rho = nothing
    l_local = nothing
    has_so = has_tag(io, "<PP_ADDINFO>")
    has_gipaw = has_tag(io, "<PP_GIPAW_RECONSTRUCTION_DATA>")
    relativistic = has_so ? "full" : "scalar"
    has_wfc = false

    return UpfHeader(generated, author, date, comment, element, pseudo_type, relativistic,
                     is_ultrasoft, is_paw, is_coulomb, has_so, has_wfc, has_gipaw,
                     paw_as_gipaw,
                     core_correction, functional, z_valence, total_psenergy, wfc_cutoff,
                     rho_cutoff, l_max, l_max_rho, l_local, mesh_size, number_of_wfc,
                     number_of_proj)
end

function upf1_parse_mesh(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_R>")
    r = read_mesh_data(Float64, io, mesh_size)

    read_until(io, "<PP_RAB>")
    rab = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)

    mesh = nothing
    rmax = nothing
    dx = nothing
    xmin = nothing
    zmesh = nothing

    return UpfMesh(r, rab, mesh, rmax, dx, xmin, zmesh)
end

function upf1_parse_nlcc(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_NLCC>")
    nlcc = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return nlcc
end

function upf1_parse_local(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_LOCAL>")
    local_ = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return local_
end

function upf1_parse_betas(io::IO, number_of_proj::Int)
    pos = position(io)
    betas = Vector{UpfBeta}(undef, number_of_proj)
    for i in 1:number_of_proj
        read_until(io, "<PP_BETA>")

        s = split(readline(io))
        index = parse(Int, s[1])
        angular_momentum = parse(Int, s[2])

        s = readline(io)
        cutoff_radius_index = parse(Int, s)

        beta = read_mesh_data(Float64, io, cutoff_radius_index)

        cutoff_radius = nothing
        norm_conserving_radius = nothing
        ultrasoft_cutoff_radius = nothing
        label = nothing

        beta = UpfBeta(beta, index, angular_momentum, cutoff_radius_index, cutoff_radius,
                       norm_conserving_radius, ultrasoft_cutoff_radius, label)
        betas[i] = beta
    end
    seek(io, pos)
    return betas
end

function upf1_parse_dij(io::IO, number_of_proj::Int)
    pos = position(io)
    read_until(io, "<PP_DIJ>")

    dij = zeros(Float64, number_of_proj, number_of_proj)

    s = split(readline(io))
    nd = parse(Int, s[1])  # Number of non-zero Dij components

    for _ in 1:nd
        s = split(readline(io))
        i, j = parse.(Int, s[1:2])
        d = parse(Float64, s[3])
        dij[i, j] = d
        dij[j, i] = d
    end
    seek(io, pos)
    return dij
end

function upf1_parse_augmentation(io::IO, mesh_size::Int, number_of_proj::Int, l_max::Int)
    pos = position(io)
    read_until(io, "<PP_QIJ>")

    s = split(readline(io))
    nqf = parse(Int, s[1])

    if nqf == 0
        return nothing
    end

    read_until(io, "<PP_RINNER>")
    rinner = Vector{Float64}(undef, 2l_max + 1)
    for i in 1:(2l_max + 1)
        s = split(readline(io))
        rinner[i] = parse(Float64, s[2])
    end
    read_until(io, "</PP_RINNER>")

    q = zeros(Float64, number_of_proj, number_of_proj)
    qijs = UpfQij[]
    qfcoefs = UpfQfcoef[]

    for i in 1:number_of_proj, _ in i:number_of_proj
        s = split(readline(io))
        first_index, second_index, second_index_angular_momentum = parse.(Int, s[1:3])

        s = split(readline(io))
        Q_int = parse(Float64, s[1])  # integral of augmentation charge

        qij = read_mesh_data(Float64, io, mesh_size)

        read_until(io, "<PP_QFCOEF>")
        qfcoef = read_mesh_data(Float64, io, nqf * (2l_max + 1))
        read_until(io, "</PP_QFCOEF>")

        is_null = nothing
        composite_index = nothing

        q[first_index, second_index] = Q_int
        q[second_index, first_index] = Q_int
        push!(qijs, UpfQij(qij, first_index, second_index, composite_index, is_null))
        push!(qfcoefs, UpfQfcoef(qfcoef, first_index, second_index, composite_index))
    end

    multipoles = nothing
    qijls = nothing
    q_with_l = false
    nqlc = nothing
    shape = nothing
    iraug = nothing
    raug = nothing
    l_max_aug = nothing
    augmentation_epsilon = nothing
    cutoff_r = nothing
    cutoff_r_index = nothing

    augmentation = UpfAugmentation(q, multipoles, qfcoefs, rinner, qijs, qijls, q_with_l,
                                   nqf, nqlc, shape, iraug, raug, l_max_aug,
                                   augmentation_epsilon, cutoff_r, cutoff_r_index)
    seek(io, pos)
    return augmentation
end

function upf1_parse_nonlocal(io::IO, pseudo_type::String, mesh_size::Int,
                             number_of_proj::Int, l_max::Int)
    betas = upf1_parse_betas(io, number_of_proj)
    dij = upf1_parse_dij(io, number_of_proj)
    if pseudo_type in ("US", "PAW")
        augmentation = upf1_parse_augmentation(io, mesh_size, number_of_proj, l_max)
    else
        augmentation = nothing
    end
    return UpfNonlocal(betas, dij, augmentation)
end

function upf1_parse_pswfc(io::IO, mesh_size::Int, number_of_wfc::Int)
    pos = position(io)
    read_until(io, "<PP_PSWFC>")
    pswfc = Vector{UpfChi}(undef, number_of_wfc)
    for index in 1:number_of_wfc
        s = split(readline(io))
        label = s[1]
        l = parse(Int, s[2])
        occupation = parse(Float64, s[3])
        
        chi = read_mesh_data(Float64, io, mesh_size)

        n = nothing
        pseudo_energy = nothing
        cutoff_radius = nothing
        ultrasoft_cutoff_radius = nothing

        pswfc[index] = UpfChi(chi, l, occupation, index, label, n, pseudo_energy, cutoff_radius,
                          ultrasoft_cutoff_radius)
    end
    seek(io, pos)
    return pswfc
end

function upf1_parse_rhoatom(io::IO, mesh_size::Int)
    pos = position(io)
    read_until(io, "<PP_RHOATOM>")
    rhoatom = read_mesh_data(Float64, io, mesh_size)
    seek(io, pos)
    return rhoatom
end

function upf1_parse_addinfo(io::IO, number_of_proj::Int, number_of_wfc::Int)
    pos = position(io)
    read_until(io, "<PP_ADDINFO>")
    relwfcs = Vector{UpfRelWfc}(undef, number_of_wfc)
    for index in 1:number_of_wfc
       s = split(readline(io))
       els = strip(s[1])
       nn = parse(Int, s[2])
       lchi = parse(Int, s[3])
       jchi = parse(Float64, s[4])
       oc = parse(Float64, s[5])
       relwfcs[index] = UpfRelWfc(jchi, index, els, nn, lchi, oc) 
    end
    relbetas = Vector{UpfRelBeta}(undef, number_of_proj)
    for index in 1:number_of_proj
        s = split(readline(io))
        lll = parse(Int, s[1])
        jjj = parse(Float64, s[2])
        relbetas[index] = UpfRelBeta(index, jjj, lll)
    end
    seek(io, pos)
    return UpfSpinOrb(relwfcs, relbetas)
end

function upf1_parse_psp(io::IO)
    version = "1.old"
    info = upf1_parse_info(io)
    header = upf1_parse_header(io)
    mesh = upf1_parse_mesh(io, header.mesh_size)
    if header.core_correction
        nlcc = upf1_parse_nlcc(io, header.mesh_size)
    else
        nlcc = nothing
    end
    if !header.is_coulomb
        local_ = upf1_parse_local(io, header.mesh_size)
    else
        local_ = nothing
    end
    nonlocal = upf1_parse_nonlocal(io, header.pseudo_type, header.mesh_size,
                                   header.number_of_proj, header.l_max)
    pswfc = upf1_parse_pswfc(io, header.mesh_size, header.number_of_wfc)
    full_wfc = nothing
    rhoatom = upf1_parse_rhoatom(io, header.mesh_size)
    if has_tag(io, "<PP_ADDINFO>")
        spinorb = upf1_parse_addinfo(io, header.number_of_proj, header.number_of_wfc)
    else
        spinorb = nothing
    end
    paw = nothing
    gipaw = nothing
    return UpfFile(version, info, header, mesh, nlcc, local_, nonlocal, pswfc, full_wfc,
                   rhoatom, spinorb, paw, gipaw)
end
