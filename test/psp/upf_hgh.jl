#TODO work on tightening tolerances
@testset "UPF-HGH agreement" begin
    integrator =
    @testset "$(splitpath(hgh_path)[end])" for (upf2_path, hgh_path) in UPF2_HGH_FILEPATHS
        upf2 = load_psp(upf2_path)
        hgh = load_psp(hgh_path)

        @testset "Local potential agrees in Fourier space" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                @test psp_quantity_evaluator(upf2, NumericLocalPotential(), FourierSpace())(q) ≈
                      psp_quantity_evaluator(hgh, NumericLocalPotential(), FourierSpace())(q) rtol = 1e-5 atol = 1e-5
            end
        end

        @testset "Nonlocal projectors agree in Fourier space" begin
            lmax = min(max_angular_momentum(upf2), max_angular_momentum(hgh))
            @testset for q in (0.0, 0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                for l in 0:lmax
                    nmax = min(n_radials(upf2, NumericProjector(), l), n_radials(hgh, NumericProjector(), l))
                    for n in 1:nmax
                        @test psp_quantity_evaluator(upf2, NumericProjector(), l, n, FourierSpace())(q) ≈
                              psp_quantity_evaluator(hgh, NumericProjector(), l, n, FourierSpace())(q) atol = 1e-6 rtol = 1e-6
                    end
                end
            end
        end

        @testset "Projector coupling coefficients agree" begin
            lmax = min(max_angular_momentum(upf2), max_angular_momentum(hgh))
            for l in 0:lmax
                nmax = min(n_radials(upf2, NumericProjector(), l), n_radials(hgh, NumericProjector(), l))
                for n in 1:nmax, m in n:nmax
                    @test get_quantity(upf2, BetaCoupling(), l, n, m) ≈
                          get_quantity(hgh, BetaCoupling(), l, n, m) atol = 1e-8 rtol = 1e-8
                end
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test psp_energy_correction(Float64, upf2) ≈
                  psp_energy_correction(Float64, hgh) rtol = 1e-5 atol = 1e-5
        end
    end
end
