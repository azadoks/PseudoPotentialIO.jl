import PseudoPotentialIO: simpson
#TODO work on tightening tolerances
@testset "UPF-HGH agreement" begin
    integrator = 
    @testset "$(splitpath(hgh_path)[end])" for (upf2_path, hgh_path) in UPF2_HGH_FILEPATHS
        upf2 = load_psp(upf2_path)
        hgh = load_psp(hgh_path)
        integrator = simpson

        @testset "Local potential agrees in Fourier space" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                @test psp_quantity_evaluator(FourierSpace(), LocalPotential(), upf2; integrator)(q) ≈
                      psp_quantity_evaluator(FourierSpace(), LocalPotential(), hgh)(q) rtol = 1e-5 atol = 1e-5
            end
        end

        @testset "Nonlocal projectors agree in Fourier space" begin
            lmax = min(max_angular_momentum(upf2), max_angular_momentum(hgh))
            @testset for q in (0.0, 0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                for l in 0:lmax
                    nmax = min(n_radials(BetaProjector(), upf2, l), n_radials(BetaProjector(), hgh, l))
                    for n in 1:nmax
                        @test psp_quantity_evaluator(FourierSpace(), BetaProjector(), upf2, l, n; integrator)(q) ≈
                              psp_quantity_evaluator(FourierSpace(), BetaProjector(), hgh, l, n)(q) atol = 1e-6 rtol = 1e-6
                    end
                end
            end
        end

        @testset "Projector coupling coefficients agree" begin
            lmax = min(max_angular_momentum(upf2), max_angular_momentum(hgh))
            for l in 0:lmax
                nmax = min(n_radials(BetaProjector(), upf2, l), n_radials(BetaProjector(), hgh, l))
                for n in 1:nmax, m in n:nmax
                    @test get_quantity(BetaCoupling(), upf2, l, n, m) ≈
                          get_quantity(BetaCoupling(), hgh, l, n, m) atol = 1e-8 rtol = 1e-8
                end
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test psp_energy_correction(Float64, upf2; integrator) ≈
                  psp_energy_correction(Float64, hgh) rtol = 1e-5 atol = 1e-5
        end
    end
end
