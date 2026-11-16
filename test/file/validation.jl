using PseudoPotentialIO: with_validation, VALIDATION, UpfHeader

# Each format: a file whose value breaks the specification is an error, a warning (and the
# value is kept), or nothing, by the validation mode
function read_with(read, text, validation)
    return read(IOBuffer(text); identifier="bad", validation)
end

function test_validation_modes(read, text, message::Regex)
    @test_throws ArgumentError read_with(read, text, :error)
    file = @test_logs (:warn, message) read_with(read, text, :warn)
    @test file isa PsPFile
    @test_logs read_with(read, text, :none)
    return file
end

@testset "Validation" begin
    @testset "HGH" begin
        text = replace(HGH_TEXTS["c-q4.hgh"], "0.34883045" => "-0.34883045")
        file = test_validation_modes(HghFile, text, r"HGH \(bad\): rloc = -0.34883045")
        @test file.rloc == -0.34883045
        # the constructor checks too
        good = load_hgh("c-q4.hgh")
        fields = [getfield(good, f) for f in fieldnames(HghFile)]
        fields[findfirst(==(:nloc), fieldnames(HghFile))] = 3
        @test_throws ArgumentError with_validation(() -> HghFile(fields...), :error)
        @test (@test_logs (:warn, r"nloc = 3") HghFile(fields...)).nloc == 3
    end

    @testset "PSP8" begin
        text = replace(PSP8_TEXTS["nlcc"], "3.0000  2.0000" => "1.0000  2.0000")
        file = test_validation_modes(Psp8File, text, r"zion = 2.0, zatom = 1.0")
        @test file.header.zatom == 1.0
        # an unknown functional: no functional, rather than an error
        text = replace(PSP8_TEXTS["nlcc"], "8   11   1" => "8   999   1")
        file = test_validation_modes(Psp8File, text, r"pspxc = 999")
        @test isempty(functional(file))
    end

    @testset "UPF" begin
        text = read(joinpath(DATA_DIR, "al.blps.lda.upf"), String)
        bad = replace(text, "relativistic=\"nonrelativistic\"" => "relativistic=\"bogus\"")
        @test bad != text
        file = test_validation_modes(UpfFile, bad, r"relativistic = \"bogus\"")
        @test file.header.relativistic == "bogus"
        # the struct is checked against the file: PP_NLCC without core_correction
        file = UpfFile(IOBuffer(text))
        h = file.header
        header = UpfHeader((f === :core_correction ? true : getfield(h, f)
                            for f in fieldnames(UpfHeader))...)
        fields = [f === :header ? header : getfield(file, f) for f in fieldnames(UpfFile)]
        @test_throws ArgumentError with_validation(() -> UpfFile(fields...), :error)
        @test with_validation(() -> UpfFile(fields...), :none).header.core_correction
    end

    # looser than the specification: a warning, never an error
    @testset "Soft checks" begin
        text = read(joinpath(DATA_DIR, "al.blps.lda.upf"), String)
        @test occursin(r"mesh_size=\"\s*\d+\"", text)
        bad = replace(text, r"mesh_size=\"\s*\d+\"" => "mesh_size=\"7\""; count=1)
        for validation in (:error, :warn)
            file = @test_logs (:warn, r"mesh_size = 7") read_with(UpfFile, bad, validation)
            @test file.header.mesh_size == 7
        end
        @test_logs read_with(UpfFile, bad, :none)
    end

    @testset "PSML" begin
        filepath = first(PSML_FILEPATHS)
        text = read(filepath, String)
        bad = replace(text, r"relativity=\"\w+\"" => "relativity=\"bogus\""; count=1)
        @test bad != text
        file = test_validation_modes(PsmlFile, bad, r"relativity = \"bogus\"")
        @test file.pseudo_atom_spec.relativity == "bogus"
    end

    @testset "Mode" begin
        text = replace(HGH_TEXTS["c-q4.hgh"], "0.34883045" => "-0.34883045")
        @test VALIDATION[] === :warn
        @test_logs (:warn, r"rloc") HghFile(IOBuffer(text))
        # scoped: for the block, and the tasks it starts
        task = with_validation(:error) do
            Threads.@spawn try HghFile(IOBuffer(text)) catch e; e end
        end
        @test fetch(task) isa ArgumentError
        # load_psp_file
        path = joinpath(mktempdir(; cleanup=true), "c.hgh")
        write(path, text)
        @test_throws ArgumentError load_psp_file(path; validation=:error)
        @test (@test_logs load_psp_file(path; validation=:none)) isa HghFile
        @test_throws ArgumentError load_psp_file(path; validation=:maybe)
        # what cannot be read is an error in any mode
        @test_throws ErrorException load_psp_file(path * "x")
        @test_throws ErrorException HghFile(IOBuffer("C\n  x"); validation=:none)
    end
end
