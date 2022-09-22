module UPF2

using EzXML

export parse_upf2

parse_bool(s::AbstractString) = occursin("T", uppercase(s)) ? true : false
parse_bool(s::Char) = uppercase(s) == 'T' ? true : false

get_attr(::Type{T}, node::EzXML.Node, key; default=zero(T)) where {T <: Number} = haskey(node, key) ? parse(T, strip(node[key])) : default
get_attr(::Type{T}, node::EzXML.Node, key; default=T("")) where {T <: AbstractString} = haskey(node, key) ? T(strip(node[key])) : default
get_attr(::Type{Bool}, node::EzXML.Node, key; default=false) = haskey(node, key) ? parse_bool(strip(node[key])) : default

function get_content(::Type{T}, node::EzXML.Node, dims...) where {T <: Number}
    text = split(strip(node.content))
    @assert length(text) == prod(dims)
    value = Array{T}(undef, dims...)
    for i = 1:prod(dims)
        value[i] = parse(T, text[i])
    end
    return value
end

function get_content(::Type{T}, node::EzXML.Node, key, dims...) where {T <: Number}
    target_node = findfirst(node, key)
    content = ismissing(target_node) ? nothing : get_content(T, target_node, dims...)
    return content
end

function parse_header!(doc_root::EzXML.Node, upf::Dict)
    header = Dict()
    # version = get_attr(String, doc_root, "version")
    node = findfirst("PP_HEADER", doc_root)
    header["format_version"] = 2
    header["generated"] = get_attr(String, node, "generated")
    header["author"] = get_attr(String, node, "author")
    header["date"] = get_attr(String, node, "date")
    header["comment"] = get_attr(String, node, "comment")
    header["element"] = get_attr(String, node, "element")
    # pseudo_type can be:
    # NC -- Norm-conserving, fully non-local only
    # SL -- Norm-conserving, non-local and semi-local
    # US -- Ultrasoft (is_ultrasoft === true)
    # 1/r -- Coulomb (is_coulomb === true)
    # PAW -- Projector-augmented wave (is_paw === true)
    header["pseudo_type"] = get_attr(String, node, "pseudo_type")
    # relativistic can be: scalar | full | nonrelativistic
    header["relativistic"] = get_attr(String, node, "relativistic")
    header["is_ultrasoft"] = get_attr(Bool, node, "is_ultrasoft")
    header["is_paw"] = get_attr(Bool, node, "is_paw")
    header["is_coulomb"] = get_attr(Bool, node, "is_coulomb")
    header["has_so"] = get_attr(Bool, node, "has_so")
    header["has_wfc"] = get_attr(Bool, node, "has_wfc")
    header["has_gipaw"] = get_attr(Bool, node, "has_gipaw")
    header["paw_as_gipaw"] = get_attr(Bool, node, "paw_as_gipaw")
    header["core_correction"] = get_attr(Bool, node, "core_correction")
    header["functional"] = get_attr(String, node, "functional")  # Max. 25 characters
    header["z_valence"] = get_attr(Float64, node, "z_valence")
    header["total_psenergy"] = get_attr(Float64, node, "total_psenergy")
    header["wfc_cutoff"] = get_attr(Float64, node, "wfc_cutoff")
    header["rho_cutoff"] = get_attr(Float64, node, "rho_cutoff")
    header["l_max"] = get_attr(Int, node, "l_max")  # 0:3
    header["l_max_rho"] = get_attr(Int, node, "l_max_rho")
    header["l_local"] = get_attr(Int, node, "l_local")  # -1 if none, 0:(l_max - 1) otherwise
    header["mesh_size"] = get_attr(Int, node, "mesh_size")
    header["number_of_wfc"] = get_attr(Int, node, "number_of_wfc")
    header["number_of_proj"] = get_attr(Int, node, "number_of_proj")

    upf["header"] = header
end

function parse_radial_grid!(doc_root::EzXML.Node, upf::Dict)
    node = findfirst("PP_MESH/PP_R", doc_root)
    upf["radial_grid"] = parse.(Float64, split(strip(node.content)))
    
    node = findfirst("PP_MESH/PP_RAB", doc_root)
    upf["radial_grid_derivative"] = parse.(Float64, split(strip(node.content)))
end

function parse_nlcc!(doc_root::EzXML.Node, upf::Dict)
    if !upf["header"]["core_correction"]
        rho = []
    else
        node = findfirst("PP_NLCC", doc_root)
        rho = parse.(Float64, split(strip(node.content)))
    end
    upf["core_charge_density"] = rho
end

function parse_local!(doc_root::EzXML.Node, upf::Dict)
    node = findfirst("PP_LOCAL", doc_root)
    upf["local_potential"] = parse.(Float64, split(strip(node.content))) ./ 2  # Ry -> Ha
end

function parse_beta_projectors!(doc_root::EzXML.Node, upf::Dict)
    beta_projectors = []
    for i = 1:upf["header"]["number_of_proj"]
        node = findfirst("PP_NONLOCAL/PP_BETA.$i", doc_root)
        beta = Dict()
        beta["radial_function"] = parse.(Float64, split(strip(node.content)))
        beta["label"] = get_attr(String, node, "label")
        beta["angular_momentum"] = get_attr(Int, node, "angular_momentum")
        beta["cutoff_radius_index"] = get_attr(Int, node, "cutoff_radius_index")
        beta["cutoff_radius"] = get_attr(Float64, node, "cutoff_radius")
        beta["index"] = get_attr(Int, node, "index")
        if upf["header"]["has_so"]
            node_so = findfirst("PP_SPIN_ORB/PP_RELBETA.$i", doc_root)
            beta["total_angular_momentum"] = get_attr(Float64, node_so, "jjj")
        end
        push!(beta_projectors, beta)
    end
    upf["beta_projectors"] = beta_projectors
end

function parse_dij!(doc_root::EzXML.Node, upf::Dict)
    node = findfirst("PP_NONLOCAL/PP_DIJ", doc_root)
    Dij = parse.(Float64, split(strip(node.content)))
    Dij = reshape(Dij, upf["header"]["number_of_proj"], upf["header"]["number_of_proj"])
    upf["D_ion"] = Dij
end

function parse_augmentation!(doc_root::EzXML.Node, upf::Dict)
    if !upf["header"]["is_ultrasoft"]
        augmentation = []
    else
        node = findfirst("PP_NONLOCAL/PP_AUGMENTATION", doc_root)
        q_with_l = get_attr(Bool, node, "q_with_l")
        if ismissing(q_with_l)
            throw(ErrorException("Parsing `q_with_l = F` is not supported."))
        end
        
        augmentation = []
        for i = 1:upf["header"]["number_of_proj"]
            li = upf["beta_projectors"][i]["angular_momentum"]
            for j = i:upf["header"]["number_of_proj"]
                lj = upf["beta_projectors"][j]["angular_momentum"]
                for l = abs(li - lj):(li + lj)
                    if (li + lj + l) % 2 == 0
                        node = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_QIJL.$i.$j.$l", doc_root)
                        Qij = Dict()
                        Qij["radial_function"] = parse.(Float64, split(strip(node.content)))
                        Qij["i"] = i
                        Qij["j"] = j
                        Qij["angular_momentum"] = get_attr(Int, node, "angular_momentum")
                        push!(augmentation, Qij)
                    end
                end
            end
        end
    end
    upf["augmentation"] = augmentation
end

function parse_paw!(doc_root::EzXML.Node, upf::Dict)
    if !(lowercase(upf["header"]["pseudo_type"]) == "paw")
        paw_data = Dict()
    else
        node = findfirst("PP_NONLOCAL/PP_AUGMENTATION", doc_root)
        upf["header"]["cutoff_radius_index"] = get_attr(Int, node, "cutoff_r_index")

        paw_data = Dict()

        node_q = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_Q", doc_root)
        paw_data["aug_integrals"] = parse.(Float64, split(strip(node_q.content)))

        node_mp = findfirst("PP_NONLOCAL/PP_AUGMENTATION/PP_MULTIPOLES", doc_root)
        paw_data["aug_multipoles"] = parse.(Float64, split(strip(node_mp.content)))

        paw_data["ae_wfc"] = []
        for i = 1:upf["header"]["number_of_proj"]
            wfc = Dict()
            node_wfc = findfirst("PP_FULL_WFC/PP_AEWFC.$i", doc_root)
            wfc["radial_function"] = parse.(Float64, split(strip(node_wfc.content)))
            wfc["angular_momentum"] = get_attr(Int, node_wfc, "l")
            wfc["label"] = get_attr(String, node_wfc, "label")
            wfc["index"] = get_attr(Int, node_wfc, "index")
            push!(paw_data["ae_wfc"], wfc)
        end

        paw_data["ps_wfc"] = []
        for i = 1:upf["header"]["number_of_proj"]
            wfc = Dict()
            node_wfc = findfirst("PP_FULL_WFC/PP_PSWFC.$i", doc_root)
            wfc["radial_function"] = parse.(Float64, split(strip(node_wfc.content)))
            wfc["angular_momentum"] = get_attr(Int, node_wfc, "l")
            # wfc["label"] = get_attr(String, node_wfc, "label")
            # wfc["index"] = get_attr(Int, node_wfc, "index")
            push!(paw_data["ps_wfc"], wfc)
        end

        node_paw = findfirst("PP_PAW", doc_root)
        paw_core_energy = get_attr(Float64, node_paw, "core_energy"; default=missing)
        if ismissing(paw_core_energy)
            # @warn "`PP_PAW` has no `core_energy` set"
        else
            upf["header"]["paw_core_energy"] = paw_core_energy ./ 2  # Ry -> Ha
        end

        node_occ = findfirst("PP_PAW/PP_OCCUPATIONS", doc_root)
        paw_data["occupations"] = parse.(Float64, split(strip(node_occ.content)))

        node_ae_nlcc = findfirst("PP_PAW/PP_AE_NLCC", doc_root)
        paw_data["ae_core_charge_density"] = parse.(Float64, split(strip(node_ae_nlcc.content)))

        node_ae_vloc = findfirst("PP_PAW/PP_AE_VLOC", doc_root)
        paw_data["ae_local_potential"] = parse.(Float64, split(strip(node_ae_vloc.content))) ./ 2  # Ry -> Ha
    end
    upf["paw_data"] = paw_data
end

function parse_pswfc!(doc_root::EzXML.Node, upf::Dict)
    atomic_wave_functions = []
    for i = 1:upf["header"]["number_of_wfc"]
        node = findfirst("PP_PSWFC/PP_CHI.$i", doc_root)
        wfc = Dict()
        wfc["label"] = get_attr(String, node, "label")
        wfc["angular_momentum"] = get_attr(Int, node, "l")
        wfc["occupation"] = get_attr(Float64, node, "occupation")
        wfc["pseudo_energy"] = get_attr(Float64, node, "pseudo_energy")
        wfc["index"] = get_attr(Int, node, "index")
        wfc["radial_function"] = parse.(Float64, split(strip(node.content)))
        if upf["header"]["has_so"]
            node_so = findfirst("PP_SPIN_ORB/PP_RELWFC.$i", doc_root)
            wfc["total_angular_momentum"] = get_attr(Float64, node_so, "jchi")
        end
        push!(atomic_wave_functions, wfc)
    end
    upf["atomic_wave_functions"] = atomic_wave_functions
end

function parse_rhoatom!(doc_root::EzXML.Node, upf::Dict)
    node = findfirst("PP_RHOATOM", doc_root)
    upf["total_charge_density"] = parse.(Float64, split(strip(node.content)))
end

function parse_spin_orbit!(doc_root::EzXML.Node, upf::Dict)
    if upf["header"]["has_so"]
        for i = 1:upf["header"]["number_of_proj"]
            node = findfirst("PP_SPIN_ORB/PP_RELBETA.$i", doc_root)
            upf["beta_projectors"][i]["angular_momentum"] = get_attr(Float64, node, "lll")
            upf["beta_projectors"][i]["total_angular_momentum"] = get_attr(Float64, node, "jjj")
        end

        # for i = 1:upf["header"]["number_of_wfc"]
        #     node = findfirst("PP_SPIN_ORB/PP_RELWFC.$i", doc_root)
        #     upf["paw_data"]["ae_wfc"]["ae_wfc_rel"] = parse(Float64, strip(node.content))
        #     upf["paw_data"]["ae_wfc"]["total_angular_momentum"] = get_attr(Float64, node, "jchi")
        #     upf["paw_data"]["ps_wfc"]["total_angular_momentum"] = get_attr(Float64, node, "jchi")
        # end
    end
end

function parse_upf2(doc_root::EzXML.Node)
    upf = Dict()

    parse_header!(doc_root, upf)
    if upf["header"]["pseudo_type"] == "PAW"
        @warn "PAW is not implemented."
    elseif upf["header"]["pseudo_type"] == "US"
        @warn "Ultrasoft is not implemented."
    end
    parse_radial_grid!(doc_root, upf)
    parse_nlcc!(doc_root, upf)
    parse_local!(doc_root, upf)
    parse_beta_projectors!(doc_root, upf)
    parse_dij!(doc_root, upf)
    # parse_augmentation!(doc_root, upf)
    # parse_paw!(doc_root, upf)
    parse_pswfc!(doc_root, upf)
    parse_rhoatom!(doc_root, upf)
    parse_spin_orbit!(doc_root, upf)

    return upf
end
end