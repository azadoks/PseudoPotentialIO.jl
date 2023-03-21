#TODO work on tightening tolerances
@testset "UPF-PSP8 agreement" begin
    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in UPF2_PSP8_FILEPATHS
        upf2 = load_psp(upf2_path)
        psp8 = load_psp(psp8_path)

        @testset "Local potential agrees in Fourier space" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0)
                @test psp_quantity_evaluator(FourierSpace(), LocalPotential(), upf2)(q) ≈
                      psp_quantity_evaluator(FourierSpace(), LocalPotential(), psp8)(q) rtol = 1e-3 atol = 1e-3
            end
        end

        @testset "Core charge density agrees in Fourier space" begin
            if !isnothing(upf2.ρcore)
                @test !isnothing(psp8.ρcore)
            end
            if !isnothing(psp8.ρcore)
                @test !isnothing(upf2.ρcore)
            end
            if !isnothing(upf2.ρcore) & !isnothing(psp8.ρcore)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0)
                    @test psp_quantity_evaluator(FourierSpace(), CoreChargeDensity(), upf2)(q) ≈
                          psp_quantity_evaluator(FourierSpace(), CoreChargeDensity(), psp8)(q) rtol = 1e-3 atol = 1e-3
                end
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test psp_energy_correction(Float64, upf2) ≈
                  psp_energy_correction(Float64, psp8) rtol = 1e-2 atol = 1e-2
        end
    end
end
