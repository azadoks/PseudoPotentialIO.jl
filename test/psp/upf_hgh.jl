@testset "UPF-HGH agreement" begin
    @testset "$(splitpath(hgh_path)[end])" for (upf2_path, hgh_path) in upf2_hgh_pairs
        upf2 = load_psp(upf2_path)
        hgh = load_psp(hgh_path)

        @testset "Local potential agrees in Fourier space" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                @test local_potential_fourier(upf2, q) ≈ local_potential_fourier(hgh, q) atol=1e-4
            end
        end

        @testset "Nonlocal projectors agree in Fourier space" begin
            @testset for q in (0.0, 0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
            for l in 0:upf2.lmax, n in eachindex(upf2.β[l])    
                    @test projector_fourier(upf2, l, n, q) ≈ projector_fourier(hgh, l, n, q) atol=1e-6
                end
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test pseudo_energy_correction(upf2) ≈ pseudo_energy_correction(hgh) rtol=1e-4 atol=1e-4
        end
    end
end