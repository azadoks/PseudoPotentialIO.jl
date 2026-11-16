#= HGH dictionaries (`file_dict.jl`) <-> `HghFile`: `_from_dict` builds the struct, with the
sizes and the symmetric coupling matrices; `_to_dict` is its inverse. =#

function _from_dict(::Type{HghFile}, d::AbstractDict; identifier="")
    # the file does not say: guess from its name
    if occursin("lda", identifier)
        functional = "LDA_C LDA_X_PZ"
    elseif occursin("pbe", identifier)
        functional = "GGA_C_PBE GGA_X_PBE"
    else
        functional = ""
    end
    h = map(d["nonlocal"]) do channel
        rows = channel["h"]
        n = length(rows)
        [rows[min(i, j)][abs(j - i) + 1] for i in 1:n, j in 1:n]
    end
    return HghFile(identifier, functional, d["title"], d["zion"], d["rloc"],
                   length(d["cloc"]), d["cloc"], length(d["nonlocal"]) - 1,
                   [channel["rp"] for channel in d["nonlocal"]], h)
end

function _to_dict(file::HghFile)
    nonlocal = [HghDict("rp" => rp, "h" => [h[i, i:end] for i in axes(h, 1)])
                for (rp, h) in zip(file.rp, file.h)]
    return HghDict("title" => file.title, "zion" => file.zion, "rloc" => file.rloc,
                   "cloc" => file.cloc, "nonlocal" => nonlocal)
end
