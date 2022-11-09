@testset "Norm conserving" begin
    @testset "From UPF" begin
        @testset "B_pbe_v1.01.uspp.F.upf" begin
            file = load_psp_file("./data/upf1/B_pbe_v1.01.uspp.F.upf")
            @test isa(load_psp(file), UltrasoftPsP)
            @test isa(load_psp("./data/upf1/B_pbe_v1.01.uspp.F.upf"), UltrasoftPsP)
        end
        @testset "Dy.GGA-PBE-paw-v1.0.UPF" begin
            file = load_psp_file("./data/upf2/Dy.GGA-PBE-paw-v1.0.UPF")
            @test_throws "Provided `UpfFile` is not an ultra" UltrasoftPsP(file)
        end
        @testset "Mg_nc-fr-04_pbesol_stringent.upf" begin
            file = load_psp_file("./data/upf2/Mg_nc-fr-04_pbesol_stringent.upf")
            @test_throws "Provided `UpfFile` is not an ultra" UltrasoftPsP(file)
        end
    end
end
