@testset "Numeric" begin
    filepaths = [
        values(psp8_filepahts)...,
        values(upf1_filepaths)...,
        upf2_filepaths["Mg.upf"],
        upf2_filepaths["He.pbe-hgh.UPF"],
        upf2_filepaths["Si.pbe-n-rrkjus_psl.1.0.0.upf"]
    ]

    @testset for filepath in filepaths
        psp = load_psp(filepath)

        @testset "Local potential Fourier agrees with real" begin
            itp = interpolate((psp.r,), psp.Vloc, (Gridded(Linear()),))
            function integral(q, r)
                return 4π * r * sphericalbesselj(0, q * r) * (r * itp(r) + psp.Zval)
            end
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                ref = quadgk(r -> integral(q, r), first(psp.r), last(psp.r))[1] -
                    4π * psp.Zval / q^2
                @test ref ≈ local_potential_fourier(psp, q) rtol = 1e-3 atol = 1e-3
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in 0:(psp.lmax), n in eachindex(psp.β[l])
                itp = interpolate((psp.r,), psp.β[l][n], (Gridded(Linear()),))
                integral(q, r) = 4π * r^2 * sphericalbesselj(l, q * r) * itp(r)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> integral(q, r), first(psp.r), last(psp.r))[1]
                    @test ref ≈ projector_fourier(psp, l, n, q) rtol = 1e-3 atol = 1e-3
                end
            end
        end

        @testset "Core charge density Fourier agrees with real" begin
            if has_nlcc(psp)
                itp = interpolate((psp.r,), psp.ρcore, (Gridded(Linear()),))
                integral(q, r) = 4π * r^2 * sphericalbesselj(0, q * r) * itp(r)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> integral(q, r), first(psp.r), last(psp.r))[1]
                    @test ref ≈ core_charge_density_fourier(psp, q) rtol = 1e-3 atol = 1e-3
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = local_potential_fourier(psp, q_small) + 4π * psp.Zval / q_small^2
            @test ref ≈ pseudo_energy_correction(psp) atol = 1e-2
        end
    end
end