# Write every UPF v2 norm-conserving and every psp8 fixture, read it back, and compare every
# field of the two structs (floating-point data to the precision of the writers).

# Paths of the fields that differ between `a` and `b`, recursively through structs and arrays
function field_mismatches(a, b; rtol=1e-12)
    mismatches = String[]
    _collect_mismatches!(mismatches, a, b, "file", rtol)
    return mismatches
end

function _collect_mismatches!(mismatches, a, b, path, rtol)
    if a isa AbstractFloat && b isa AbstractFloat
        isapprox(a, b; rtol) || push!(mismatches, path)
    elseif a isa AbstractArray{<:AbstractFloat} && b isa AbstractArray{<:AbstractFloat}
        size(a) == size(b) && isapprox(a, b; rtol) || push!(mismatches, path)
    elseif a isa AbstractDict && b isa AbstractDict
        if keys(a) != keys(b)
            push!(mismatches, path)
        else
            for k in keys(a)
                _collect_mismatches!(mismatches, a[k], b[k], "$path[$(repr(k))]", rtol)
            end
        end
    elseif a isa AbstractArray && b isa AbstractArray
        if size(a) != size(b)
            push!(mismatches, path)
        else
            for i in eachindex(a, b)
                _collect_mismatches!(mismatches, a[i], b[i], "$path[$i]", rtol)
            end
        end
    elseif isstructtype(typeof(a)) && typeof(a) == typeof(b) && !(a isa AbstractString) &&
           fieldcount(typeof(a)) > 0
        for n in fieldnames(typeof(a))
            _collect_mismatches!(mismatches, getfield(a, n), getfield(b, n), "$path.$n", rtol)
        end
    else
        isequal(a, b) || push!(mismatches, path)
    end
    return mismatches
end

function roundtrip(file, ext)
    outpath = joinpath(mktempdir(; cleanup=true), "roundtrip" * ext)
    save_psp_file(outpath, file)
    return load_psp_file(outpath)
end

# the identifier is the file name
without_identifier(file::UpfFile) = [getfield(file, n) for n in fieldnames(UpfFile)
                                     if n != :identifier]
without_identifier(file::Psp8File) = [getfield(file, n) for n in fieldnames(Psp8File)
                                      if n != :identifier]

@testset "Round trip" begin
    # UPF v1 and v2 ultrasoft and PAW, written as UPF v2: augmentation (PP_QIJ, PP_QIJL,
    # PP_QFCOEF, PP_RINNER) and the PP_AEWFC/PP_PSWFC of PP_FULL_WFC and PP_PAW
    @testset "UPF US/PAW: $name" for (name, filepath) in
                                     merge(UPF1_CASE_FILEPATHS, UPF2_CASE_FILEPATHS)
        file = load_psp_file(filepath)
        is_norm_conserving(file) && continue
        newfile = roundtrip(file, ".upf")
        compared = [n for n in fieldnames(UpfFile) if !(n in (:identifier, :version))]
        mismatches = field_mismatches([getfield(file, n) for n in compared],
                                      [getfield(newfile, n) for n in compared])
        @test isempty(mismatches) || (filepath, mismatches)
    end

    @testset "UPF v2.0.1 NC" begin
        for filepath in values(UPF2_CASE_FILEPATHS)
            file = load_psp_file(filepath)
            is_norm_conserving(file) || continue
            newfile = roundtrip(file, ".upf")
            mismatches = field_mismatches(without_identifier(file),
                                          without_identifier(newfile))
            @test isempty(mismatches) || (filepath, mismatches)
        end
    end

    # the radii of the kinks of the tables (an extension of UPF): written as attributes of
    # the elements, with the number of mesh points below each, and read back
    @testset "UPF breakpoints" begin
        filepath = first(p for p in values(UPF2_CASE_FILEPATHS)
                         if is_norm_conserving(load_psp_file(p)) &&
                            !isnothing(load_psp_file(p).local_) &&
                            !isempty(load_psp_file(p).nonlocal.betas) &&
                            !isempty(load_psp_file(p).pswfc))
        file = load_psp_file(filepath)
        r = file.mesh.r
        breaks = Dict("PP_LOCAL" => [r[200] + 1e-3, r[300] / 3],
                      "PP_BETA.$(first(file.nonlocal.betas).index)" => [r[150] + 1e-6],
                      "PP_RHOATOM" => [r[120] + 1e-4],
                      "PP_CHI.$(first(file.pswfc).index)" => [r[180] + 2e-5])
        fields = [getfield(file, n) for n in fieldnames(UpfFile)]
        marked = UpfFile(fields[1:(end - 1)]..., breaks)
        outpath = joinpath(mktempdir(; cleanup=true), "marked.upf")
        save_psp_file(outpath, marked)
        @test load_psp_file(outpath).breakpoints == breaks
        text = read(outpath, String)
        indices = join((searchsortedfirst(r, b) - 1 for b in breaks["PP_LOCAL"]), ' ')
        @test occursin("breakpoint_indices=\"$indices\"", text)
        @test isempty(file.breakpoints)
    end

    @testset "PSP8" begin
        for filepath in values(PSP8_FILEPATHS)
            file = load_psp_file(filepath)
            newfile = roundtrip(file, ".psp8")
            mismatches = field_mismatches(without_identifier(file),
                                          without_identifier(newfile))
            @test isempty(mismatches) || (filepath, mismatches)
        end
    end
end

@testset "Libxc functional strings" begin
    qe_to_libxc = PseudoPotentialIO.qe_libxc_to_libxc
    @test qe_to_libxc("XC-001L-012L") == "LDA_X LDA_C_PW"
    @test qe_to_libxc("XC-000I-000I-101L-130L") == "GGA_X_PBE GGA_C_PBE"
    @test qe_to_libxc("XC-000I-000I-000I-000I-497L-498L") == "MGGA_X_R2SCAN MGGA_C_R2SCAN"
    @test qe_to_libxc("XC-001I-001I") == "LDA_X LDA_C_PZ"  # QE-internal sla+pz
    @test qe_to_libxc("XC-000I-000I-000I-000I-001I") == "UNKNOWN"  # QE-internal meta-GGA
    for s in ("LDA_X LDA_C_PW", "GGA_X_PBE GGA_C_PBE", "MGGA_X_R2SCAN MGGA_C_R2SCAN")
        @test qe_to_libxc(PseudoPotentialIO.libxc_to_qe_libxc(s)) == s
    end

    @test PseudoPotentialIO.libxc_to_qe_libxc("LDA_X LDA_C_PW") == "XC-001L-012L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("GGA_X_PBE GGA_C_PBE") ==
          "XC-000I-000I-101L-130L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("GGA_C_PBE GGA_X_PBE") ==
          "XC-000I-000I-101L-130L"
    @test PseudoPotentialIO.libxc_to_qe_libxc("MGGA_X_R2SCAN MGGA_C_R2SCAN") ==
          "XC-000I-000I-000I-000I-497L-498L"
end
