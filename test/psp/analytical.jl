@testset "Analytical" begin
    function local_potential_integrand(psp::AnalyticalPsP, q)
        Vloc = local_potential_real(psp)
        return r -> 4π * r * sphericalbesselj(0, q * r) * (r * Vloc(r) + valence_charge(psp))
    end

    function projector_integrand(psp::AnalyticalPsP, l::Int, n::Int, q)
        β = projector_real(psp, l, n)
        return r -> 4π * r^2 * sphericalbesselj(l, q * r) * β(r)
    end

    rtest = 200:-0.1:0
    itest = randperm(length(HGH_FILEPATHS))[1:30]

    @testset "$(splitpath(filepath)[end])" for filepath in HGH_FILEPATHS[itest]
        psp = load_psp(filepath)

        @testset "Local potential Fourier agrees with real" begin
            Vloc = local_potential_real(psp)
            rcut = rtest[findfirst(r -> !isapprox(0, Vloc(r)), rtest)]
            Ṽloc = local_potential_fourier(psp)
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                ref = quadgk(local_potential_integrand(psp, q), 0, rcut)[1] -
                      4π * valence_charge(psp) / q^2
                @test ref ≈ Ṽloc(q) rtol = 1e-9 atol = 1e-9
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in angular_momenta(psp), n in projector_radial_indices(psp, l)
                β = projector_real(psp, l, n)
                rcut = rtest[findfirst(r -> !isapprox(0, β(r)), rtest)]
                β̃ = projector_fourier(psp, l, n)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(projector_integrand(psp, l, n, q), 0, rcut)[1]
                    @test ref ≈ β̃(q) rtol = 1e-9 atol = 1e-9
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = local_potential_fourier(psp)(q_small) +
                  4π * valence_charge(psp) / q_small^2
            @test ref ≈ pseudo_energy_correction(Float64, psp) atol = 1e-2
        end
    end
end
