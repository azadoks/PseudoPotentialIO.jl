@testset "Norm conserving" begin
    @testset "From PSP8" begin
        @testset "H_nc-sr-04_pbe_standard.psp8" begin
            psp = load_psp("./data/psp8/H_nc-sr-04_pbe_standard.psp8")
            @test isa(psp, NormConservingPsP)
        end
        @testset "Ti_nc-fr-04_pbesol_standard.psp8" begin
            @test_throws "Fully relativistic" load_psp("./data/psp8/Ti_nc-fr-04_pbesol_standard.psp8")
        end
        @testset "Zn_nc-sr-04_pbesol_stringent.psp8" begin
            psp = load_psp("./data/psp8/Zn_nc-sr-04_pbesol_stringent.psp8")
            @test isa(psp, NormConservingPsP)
        end
    end

    @testset "From UPF" begin
        @testset "B_pbe_v1.01.uspp.F.upf" begin
            file = load_psp_file("./data/upf1/B_pbe_v1.01.uspp.F.upf")
            @test_throws "Provided `UpfFile` is not a norm" NormConservingPsP(file)
        end
        @testset "Dy.GGA-PBE-paw-v1.0.UPF" begin
            file = load_psp_file("./data/upf2/Dy.GGA-PBE-paw-v1.0.UPF")
            @test_throws "Provided `UpfFile` is not a norm" NormConservingPsP(file)
        end
        @testset "Mg_nc-fr-04_pbesol_stringent.upf" begin
            @test_throws "Fully relativistic" load_psp("./data/upf2/Mg_nc-fr-04_pbesol_stringent.upf")
        end
    end
end
