@testset "UPF-PSP8 agreement" begin
    QUANTITIES = [NumericLocalPotential(), CoreChargeDensity(), ValenceChargeDensity()]
    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in
                                                UPF2_PSP8_FILEPATHS
        qs = [0.01, 0.5, 2.5, 5.0, 10.0]

        upf2 = load_psp(upf2_path)
        psp8 = load_psp(psp8_path)

        # You _must_ truncate the UPF at the same grid point as in the psp8 to get
        # results that agree. Often, it seems that this trunction is _before_
        # the quantity actually decays to zero!
        upf2_interpolated = interpolate_evaluate(upf2, psp8.r)

        upf2_q = hankel_transform(upf2_interpolated, qs)
        psp8_q = hankel_transform(psp8_interpolated, qs)

        @testset "$quantity agrees in Fourier space" for quantity in QUANTITIES
            if has_quantity(upf2, quantity) && has_quantity(psp8, quantity)
                upf2_quantity = get_quantity(upf2_q, quantity)
                psp8_quantity = get_quantity(psp8_q, quantity)

                @testset "q=$(q)" for (i, q) in enumerate(qs)
                    @test upf2_quantity[i] ≈ psp8_quantity[i] rtol = 1e-8 atol = 1e-8
                end
            end
        end

        @testset "$quantity agrees in Real space" for quantity in QUANTITIES
            if has_quantity(upf2, quantity) && has_quantity(psp8, quantity)

                upf2_quantity = get_quantity(upf2_interpolated, quantity)
                psp8_quantity = get_quantity(psp8, quantity)

                @test all(isapprox.(upf2_evaluator.(rgrid), psp8_evaluator.(rgrid),
                                    rtol=1e-9, atol=1e-9))
            end
        end

        @testset "Pseudo energy correction agrees" begin
            @test psp_energy_correction(Float64, upf2, psp8.r, psp8.dr, psp8.Δr) ≈
                  psp_energy_correction(Float64, psp8, psp8.r, psp8.dr, psp8.Δr) rtol = 1e-9 atol = 1e-9
        end
    end
end
