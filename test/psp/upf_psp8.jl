#TODO work on tightening tolerances
@testset "UPF-PSP8 agreement" begin
    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in upf2_psp8_pairs
        upf2 = load_psp(upf2_path)
        psp8 = load_psp(psp8_path)

        @testset "Local potential agrees in Fourier space" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                @test local_potential_fourier(upf2, q) ≈
                      local_potential_fourier(psp8, q) rtol = 1e-2 atol = 1e-2
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
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    @test local_potential_fourier(upf2, q) ≈
                          local_potential_fourier(psp8, q) rtol = 1e-2 atol = 1e-2
                end
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test pseudo_energy_correction(upf2) ≈
                  pseudo_energy_correction(psp8) rtol = 1e-2 atol = 1e-1
        end
    end
end
