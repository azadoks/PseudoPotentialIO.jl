@testset "UPF-PSP8 agreement" begin
    QUANTITIES = [LocalPotential(), CoreChargeDensity(), ValenceChargeDensity()]
    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in
                                                UPF2_PSP8_FILEPATHS
        upf2 = load_psp(upf2_path)
        psp8 = load_psp(psp8_path)

        @testset "$quantity agrees in Fourier space" for quantity in QUANTITIES
            if has_quantity(upf2, quantity) && has_quantity(psp8, quantity)
                # You _must_ truncate the UPF at the same grid point as in the psp8 to get
                # results that agree. Often, it seems that this trunction is _before_
                # the quantity actually decays to zero!
                upf2_evaluator = psp_quantity_evaluator(upf2, quantity, FourierSpace(), psp8.r, psp8.dr, psp8.Δr)
                psp8_evaluator = psp_quantity_evaluator(psp8, quantity, FourierSpace(), psp8.r, psp8.dr, psp8.Δr)

                @testset "q=$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0)
                    @test upf2_evaluator(q) ≈ psp8_evaluator(q) rtol = 1e-8 atol = 1e-8
                end
            end
        end

        @testset "$quantity agrees in Real space" for quantity in QUANTITIES
            if has_quantity(upf2, quantity) && has_quantity(psp8, quantity)
                rmin = max(first(upf2.r), first(psp8.r))
                rmax = min(cutoff_radius(upf2, quantity), cutoff_radius(psp8, quantity))
                rgrid = rmin:0.01:rmax

                upf2_evaluator = psp_quantity_evaluator(upf2, quantity, RealSpace())
                psp8_evaluator = psp_quantity_evaluator(psp8, quantity, RealSpace())

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
