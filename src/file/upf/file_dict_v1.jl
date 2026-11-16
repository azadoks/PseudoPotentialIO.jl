#= UPF v1 -> nested dictionaries

UPF v1 is read into the dictionaries of UPF v2 (`file_dict_v2.jl`), as QE reads both into
one `pseudo_upf`: each value under the name of its UPF v2 element or attribute (line 6 of
`PP_HEADER` is `header["z_valence"]`, `PP_ADDINFO` is "spin_orb", ...), and "version" is
"1.old". There is no writer: UPF v1 is written as UPF v2.

Based on QE's reader (`read_upf_v1.f90` in upflib):
  - a block starts at the first line whose first word contains its tag (`<PP_MESH>`), from
    the start of the file (from the current line inside `PP_MESH` and `PP_NONLOCAL`)
  - values are read list-directed: the values needed, the rest of the line is a comment;
    an array runs over as many lines as it needs, and the rest of its last line is skipped
  - the functional is the first 20 characters of its line
  - the labels, l and occupations of the wavefunctions are those of `PP_HEADER` (the label
    line of each wavefunction in `PP_PSWFC` is skipped)
  - a projector can have a line with its cutoff radii and one with its label after its
    values; it has the values up to its cutoff index, the rest are zeros
  - the augmentation functions are given for each l (`PP_QIJL`) if the file has a
    `<PP_QIJ_WITH_L>` or is PAW; the l in the first line of each is checked
  - `PP_ADDINFO` (spin-orbit, and the parameters of the mesh) is ignored if a j is not
    l ± 1/2, as in obsolete files
GIPAW data (`PP_GIPAW_RECONSTRUCTION_DATA`) is not read.
=#

"""
Read a UPF v1 file (its text) into the nested dictionaries of UPF v2 (see
`src/file/upf/file_dict_v1.jl`).
"""
function upf1_to_dict(text::AbstractString)::UpfDict
    r = FortranRecords(text)
    try
        return upf1_records_to_dict(r)
    catch e
        e isa ErrorException ? error("UPF v1, line $(r.line): $(e.msg)") : rethrow()
    end
end

upf1_first_word(line) = (words = split(line; limit=2); isempty(words) ? "" : words[1])
upf1_has(r, tag) = any(line -> occursin("<PP_$tag>", upf1_first_word(line)), r.lines)

# Go to the line of `<PP_$tag>`: the first one after the current line, or in the file
function upf1_scan!(r, tag; rewind=false)
    for i in (rewind ? 1 : r.line + 1):length(r.lines)
        if occursin("<PP_$tag>", upf1_first_word(r.lines[i]))
            r.line = i
            return
        end
    end
    error("no <PP_$tag>")
end

# `n` reals, from as many lines as they take
function upf1_reals!(r, n)
    values = Float64[]
    while length(values) < n
        for word in split(next_line!(r))
            length(values) == n && break
            x = parse_fortran_reals(word)
            isnothing(x) && error("$(repr(word)) is not a number")
            append!(values, x)
        end
    end
    length(values) == n || error("$(length(values)) values, but $n are needed")
    return values
end

upf1_real!(r) = fortran_real(only(next_words!(r, 1)))
upf1_int!(r) = fortran_int(only(next_words!(r, 1)))

# a list-directed logical: `T`, `.true.`, `F`, ...
function upf1_logical(word)
    m = match(r"^\.?([TtFf])", word)
    isnothing(m) && error("$(repr(word)) is not a logical")
    return uppercase(m[1]) == "T"
end

function upf1_records_to_dict(r::FortranRecords)
    d = UpfDict("version" => "1.old")
    if upf1_has(r, "INFO")
        upf1_scan!(r, "INFO"; rewind=true)
        stop = something(findnext(l -> occursin("</PP_INFO>", l), r.lines, r.line), r.line)
        d["info"] = UpfDict("content" => join(r.lines[(r.line + 1):(stop - 1)], '\n'))
    end
    d["header"] = header = upf1_header!(r)
    mesh = header["mesh_size"]
    wavefunctions = pop!(header, "wavefunctions")

    upf1_scan!(r, "MESH"; rewind=true)
    upf1_scan!(r, "R")
    radii = upf1_reals!(r, mesh)
    upf1_scan!(r, "RAB")
    d["mesh"] = UpfDict("mesh" => mesh, "r" => radii, "rab" => upf1_reals!(r, mesh))
    if header["core_correction"]
        upf1_scan!(r, "NLCC"; rewind=true)
        d["nlcc"] = upf1_reals!(r, mesh)
    end
    if !header["is_coulomb"]
        upf1_scan!(r, "LOCAL"; rewind=true)
        d["local"] = upf1_reals!(r, mesh)
        upf1_scan!(r, "NONLOCAL"; rewind=true)
        d["nonlocal"] = upf1_nonlocal!(r, header, upf1_has(r, "QIJ_WITH_L"))
    end
    upf1_scan!(r, "PSWFC"; rewind=true)
    chis = map(enumerate(wavefunctions)) do (index, (label, l, occupation))
        next_line!(r)  # the label, l and occupation again
        UpfDict("label" => label, "l" => l, "occupation" => occupation, "index" => index,
                "content" => upf1_reals!(r, mesh))
    end
    d["pswfc"] = UpfDict("chi" => chis)
    upf1_scan!(r, "RHOATOM"; rewind=true)
    d["rhoatom"] = upf1_reals!(r, mesh)
    if upf1_has(r, "ADDINFO")
        spin_orb, mesh_parameters = upf1_addinfo!(r, header)
        if !isnothing(spin_orb)
            header["has_so"] = true
            header["relativistic"] = "full"
            d["spin_orb"] = spin_orb
            merge!(d["mesh"], mesh_parameters)
        end
    end
    return d
end

function upf1_header!(r)
    upf1_scan!(r, "HEADER"; rewind=true)
    next_words!(r, 1)  # the version of the pseudopotential
    header = UpfDict("element" => String(only(next_words!(r, 1))))
    pseudo_type = String(only(next_words!(r, 1)))
    header["pseudo_type"] = pseudo_type
    header["relativistic"] = "scalar"
    if occursin("US", pseudo_type)
        kind = (false, true, false)
    elseif occursin("PAW", pseudo_type)
        kind = (true, true, false)
    elseif occursin("NC", pseudo_type)
        kind = (false, false, false)
    elseif occursin("1/r", pseudo_type)
        kind = (false, false, true)
    else
        error("unknown type of pseudopotential $(repr(pseudo_type))")
    end
    header["is_paw"], header["is_ultrasoft"], header["is_coulomb"] = kind
    header["has_so"] = header["has_wfc"] = header["has_gipaw"] = false
    header["core_correction"] = upf1_logical(only(next_words!(r, 1)))
    header["with_metagga_info"] = false
    header["functional"] = String(strip(first(next_line!(r), 20)))
    header["z_valence"] = upf1_real!(r)
    header["total_psenergy"] = upf1_real!(r)
    header["wfc_cutoff"], header["rho_cutoff"] = fortran_real.(next_words!(r, 2))
    header["l_max"] = upf1_int!(r)
    header["mesh_size"] = upf1_int!(r)
    header["number_of_wfc"], header["number_of_proj"] = fortran_int.(next_words!(r, 2))
    next_line!(r)  # Wavefunctions  nl  l  occ
    header["wavefunctions"] = map(1:header["number_of_wfc"]) do _
        label, l, occupation = next_words!(r, 3)
        (String(label), fortran_int(l), fortran_real(occupation))
    end
    return header
end

function upf1_nonlocal!(r, header, q_with_l)
    nbeta, mesh = header["number_of_proj"], header["mesh_size"]
    nbeta == 0 && return UpfDict()
    betas = map(1:nbeta) do index
        upf1_scan!(r, "BETA")
        l = fortran_int(next_words!(r, 2)[2])
        ikk = upf1_int!(r)
        ikk <= mesh || error("PP_BETA $index has a cutoff index $ikk beyond the mesh ($mesh)")
        beta = UpfDict("index" => index, "angular_momentum" => l, "cutoff_radius_index" => ikk)
        beta["content"] = [upf1_reals!(r, ikk); zeros(mesh - ikk)]
        # optional: the cutoff radii, then the label
        ended() = r.line == length(r.lines) || occursin("</PP_BETA>", r.lines[r.line + 1])
        if !ended() && length(split(r.lines[r.line + 1])) >= 2
            radii = parse_fortran_reals.(split(r.lines[r.line + 1])[1:2])
            if all(x -> !isnothing(x) && length(x) == 1, radii)
                next_line!(r)
                beta["cutoff_radius"], beta["ultrasoft_cutoff_radius"] = only.(radii)
                ended() || (beta["label"] = String(only(next_words!(r, 1))))
            end
        end
        return beta
    end
    upf1_scan!(r, "DIJ")
    dij = zeros(nbeta, nbeta)
    for _ in 1:upf1_int!(r)
        i, j, x = next_words!(r, 3)
        dij[fortran_int(i), fortran_int(j)] = dij[fortran_int(j), fortran_int(i)] = fortran_real(x)
    end
    d = UpfDict("beta" => betas, "dij" => vec(dij))
    header["is_ultrasoft"] &&
        (d["augmentation"] = upf1_augmentation!(r, header, betas, q_with_l || header["is_paw"]))
    return d
end

function upf1_augmentation!(r, header, betas, with_l)
    nbeta, mesh = length(betas), header["mesh_size"]
    upf1_scan!(r, "QIJ")
    nqf = upf1_int!(r)
    nqlc = 2 * header["l_max"] + 1
    d = UpfDict("q_with_l" => with_l, "nqf" => nqf, "nqlc" => nqlc)
    if nqf > 0
        upf1_scan!(r, "RINNER")
        d["rinner"] = upf1_reals!(r, 2nqlc)[2:2:end]  # index, radius
        next_line!(r)  # </PP_RINNER>
    end
    q = zeros(nbeta, nbeta)
    qfcoef = zeros(nqf, nqlc, nbeta, nbeta)
    functions = UpfDict[]
    for i in 1:nbeta, j in i:nbeta
        lj = betas[j]["angular_momentum"]
        l = fortran_int(next_words!(r, 3)[3])
        l == lj || error("Q($i, $j) has l = $l, but projector $j has l = $lj")
        q[i, j] = q[j, i] = upf1_real!(r)
        li = betas[i]["angular_momentum"]
        for L in (with_l ? (abs(li - lj):(li + lj)) : (nothing,))
            f = UpfDict("first_index" => i, "second_index" => j)
            isnothing(L) || (f["angular_momentum"] = L)
            f["content"] = upf1_reals!(r, mesh)
            push!(functions, f)
        end
        if nqf > 0
            upf1_scan!(r, "QFCOEF")
            qfcoef[:, :, i, j] = qfcoef[:, :, j, i] = reshape(upf1_reals!(r, nqf * nqlc), nqf, nqlc)
            next_line!(r)  # </PP_QFCOEF>
        end
    end
    d["q"] = vec(q)
    nqf > 0 && (d["qfcoef"] = vec(qfcoef))
    d[with_l ? "qijl" : "qij"] = functions
    return d
end

# "spin_orb" and the parameters of the mesh, or `nothing` (obsolete PP_ADDINFO)
function upf1_addinfo!(r, header)
    upf1_scan!(r, "ADDINFO"; rewind=true)
    half(j, l) = abs(abs(j - l) - 0.5) <= 1e-7
    relwfcs = map(1:header["number_of_wfc"]) do index
        els, nn, lchi, jchi, oc = next_words!(r, 5)
        UpfDict("index" => index, "els" => String(els), "nn" => fortran_int(nn),
                "lchi" => fortran_int(lchi), "jchi" => fortran_real(jchi),
                "oc" => fortran_real(oc))
    end
    all(w -> half(w["jchi"], w["lchi"]), relwfcs) || return nothing, nothing
    relbetas = map(1:header["number_of_proj"]) do index
        lll, jjj = next_words!(r, 2)
        UpfDict("index" => index, "lll" => fortran_int(lll), "jjj" => fortran_real(jjj))
    end
    all(b -> half(b["jjj"], b["lll"]), relbetas) || return nothing, nothing
    xmin, rmax, zmesh, dx = fortran_real.(next_words!(r, 4))
    return (UpfDict("relwfc" => relwfcs, "relbeta" => relbetas),
            UpfDict("xmin" => xmin, "rmax" => rmax, "zmesh" => zmesh, "dx" => dx))
end
