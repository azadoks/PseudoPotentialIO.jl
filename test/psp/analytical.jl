@testset "Analytical" begin
    function local_potential_integrand(psp::AnalyticalPsP, q)
        Vloc = psp_quantity_evaluator(psp, NumericLocalPotential(), RealSpace())
        return r -> 4π * r * sphericalbesselj(0, q * r) *
                    (r * Vloc(r) + valence_charge(psp))
    end

    function beta_projector_integrand(psp::AnalyticalPsP, l::Int, n::Int, q)
        β = psp_quantity_evaluator(psp, NumericProjector(), l, n, RealSpace())
        return r -> 4π * r^2 * sphericalbesselj(l, q * r) * β(r)
    end

    rtest = 200:-0.1:0
    itest = randperm(length(HGH_FILEPATHS))[1:30]  # ≈3s/pseudo, 237 HGH in total

    @testset "$(splitpath(filepath)[end])" for filepath in HGH_FILEPATHS[itest]
        psp = load_psp(filepath)

        @testset "Local potential Fourier agrees with real" begin
            Vloc = psp_quantity_evaluator(psp, NumericLocalPotential(), RealSpace())
            rcut = rtest[findfirst(r -> !isapprox(0, Vloc(r)), rtest)]
            Ṽloc = psp_quantity_evaluator(psp, NumericLocalPotential(), FourierSpace())
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                ref = quadgk(local_potential_integrand(psp, q), 0, rcut)[1] -
                      4π * valence_charge(psp) / q^2
                @test ref ≈ Ṽloc(q) rtol = 1e-9 atol = 1e-9
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in angular_momenta(psp), n in 1:n_radials(psp, NumericProjector(), l)
                β = psp_quantity_evaluator(psp, NumericProjector(), l, n, RealSpace())
                rcut = rtest[findfirst(r -> !isapprox(0, β(r)), rtest)]
                β̃ = psp_quantity_evaluator(psp, NumericProjector(), l, n, FourierSpace())
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(beta_projector_integrand(psp, l, n, q), 0, rcut)[1]
                    @test ref ≈ β̃(q) rtol = 1e-9 atol = 1e-9
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = psp_quantity_evaluator(psp, NumericLocalPotential(), FourierSpace())(q_small) +
                  4π * valence_charge(psp) / q_small^2
            @test ref ≈ psp_energy_correction(Float64, psp) atol = 1e-2
        end

        if max_angular_momentum(psp) >= 0
            β = psp_quantity_evaluator(psp, NumericProjector(), 0, 4, FourierSpace())
            @test_throws ArgumentError β(0.1)
        end
    end
end
