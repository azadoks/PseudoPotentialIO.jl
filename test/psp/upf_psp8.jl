@testset "UPF-PSP8 agreement" begin
    QUANTITIES = [LocalPotential(), CoreChargeDensity(), ValenceChargeDensity()]
    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in
                                                UPF2_PSP8_FILEPATHS
        upf2 = load_psp(upf2_path)
        psp8 = load_psp(psp8_path)

        @testset "$quantity agrees in Fourier space" for quantity in QUANTITIES
            if has_quantity(quantity, upf2) && has_quantity(quantity, psp8)
                # You _must_ truncate the UPF at the same grid point as in the psp8 to get
                # results that agree. Often, it seems that this trunction is _before_
                # the quantity actually decays to zero!
                i_stop = lastindex(get_quantity(quantity, psp8))
                upf2_evaluator = psp_quantity_evaluator(FourierSpace(), quantity, upf2,
                                                        i_stop)
                psp8_evaluator = psp_quantity_evaluator(FourierSpace(), quantity, psp8)

                @testset "q=$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0)
                    @test upf2_evaluator(q) ≈ psp8_evaluator(q) rtol = 1e-8 atol = 1e-8
                end
            end
        end

        @testset "$quantity agrees in Real space" for quantity in QUANTITIES
            if has_quantity(quantity, upf2) && has_quantity(quantity, psp8)
                rmin = max(first(upf2.r), first(psp8.r))
                rmax = min(cutoff_radius(quantity, upf2), cutoff_radius(quantity, psp8))
                rgrid = rmin:0.01:rmax

                upf2_evaluator = psp_quantity_evaluator(RealSpace(), quantity, upf2)
                psp8_evaluator = psp_quantity_evaluator(RealSpace(), quantity, psp8)

                @test all(isapprox.(upf2_evaluator.(rgrid), psp8_evaluator.(rgrid),
                                    rtol=1e-9, atol=1e-9))
            end
        end

        @testset "Pseudo energy correction agrees" begin
            i_stop = lastindex(get_quantity(LocalPotential(), psp8))
            @test psp_energy_correction(Float64, upf2, i_stop) ≈
                  psp_energy_correction(Float64, psp8) rtol = 1e-9 atol = 1e-9
        end
    end
end
