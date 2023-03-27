import PseudoPotentialIO: guess_mesh_type

#TODO work on tightening tolerances
#TODO * fixed trapezoid integration, allowed for tighter tolerances

#* Note: q-values are equivalent to the largest q at an Ecut via q^2 = 2Ecut
#*       Ecut up to ≈250 Ha is considered.
@testset "Numeric" begin
    function local_potential_integrand(psp::NumericPsP, q)
        Vloc = psp_quantity_evaluator(RealSpace(), LocalPotential(), psp)
        return r -> 4π * r * sphericalbesselj(0, q * r) * (r * Vloc(r) + psp.Zval)
    end

    function radial_quantity_integrand(::NumericPsP, f, l::Int, q)
        return r -> 4π * sphericalbesselj(l, q * r) * f(r)
    end

    @testset "$(name)" for (name, filepath) in NUMERIC_CASE_FILEPATHS
        psp = load_psp(filepath)
        rmin = max(first(psp.r), eps(eltype(psp.r)))

        @testset "Local potential Fourier agrees with real" begin
            rmax = cutoff_radius(LocalPotential(), psp)
            Vloc_fourier = psp_quantity_evaluator(FourierSpace(), LocalPotential(), psp)
            @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                ref = quadgk(local_potential_integrand(psp, q), rmin, rmax)[1] -
                    4π * psp.Zval / q^2
                @test ref ≈ Vloc_fourier(q) rtol = 1e-2 atol = 1e-2
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in angular_momenta(psp), n in 1:n_radials(BetaProjector(), psp, l)
                rmax = cutoff_radius(BetaProjector(), psp, l, n)
                β_real = psp_quantity_evaluator(RealSpace(), BetaProjector(), psp, l, n)
                β_fourier = psp_quantity_evaluator(FourierSpace(), BetaProjector(), psp, l, n)
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, β_real, l, q), rmin, rmax)[1]
                    @test ref ≈ β_fourier(q) rtol = 1e-2 atol = 1e-2
                end
            end
        end

        @testset "Core charge density Fourier agrees with real" begin
            if has_quantity(CoreChargeDensity(), psp)
                rmax = cutoff_radius(CoreChargeDensity(), psp)
                ρcore_real = psp_quantity_evaluator(RealSpace(), CoreChargeDensity(), psp)
                ρcore_fourier = psp_quantity_evaluator(FourierSpace(), CoreChargeDensity(), psp)
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, ρcore_real, 0, q), rmin, rmax)[1]
                    @test ref ≈ ρcore_fourier(q) rtol = 1e-4 atol = 1e-4
                end
            end
        end

        @testset "Pseudo-atomic wavefunction Fouriers agree with real" begin
            if has_quantity(ChiProjector(), psp)
                for l in angular_momenta(psp), n in 1:n_radials(ChiProjector(), psp, l)
                    rmax = cutoff_radius(ChiProjector(), psp, l, n)
                    χ_real = psp_quantity_evaluator(RealSpace(), ChiProjector(), psp, l, n)
                    χ_fourier = psp_quantity_evaluator(FourierSpace(), ChiProjector(), psp, l, n)
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, χ_real, l, q), rmin, rmax)[1]
                        @test ref ≈ χ_fourier(q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end
        end

        @testset "Valence charge density Fouriers agree with real" begin
            if has_quantity(ValenceChargeDensity(), psp)
                rmax = cutoff_radius(ValenceChargeDensity(), psp)
                ρval_real = psp_quantity_evaluator(RealSpace(), ValenceChargeDensity(), psp)
                ρval_fourier = psp_quantity_evaluator(FourierSpace(), ValenceChargeDensity(), psp)
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, ρval_real, 0, q), rmin, rmax)[1]
                    @test ref ≈ ρval_fourier(q) rtol = 1e-4 atol = 1e-4
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = psp_quantity_evaluator(FourierSpace(), LocalPotential(), psp)(q_small) + 4π * psp.Zval / q_small^2
            @test ref ≈ psp_energy_correction(Float64, psp) rtol = 1e-4 atol = 1e-4
        end

        # Test resampling log-mesh pseudopotentials on a linear mesh using cubic splines before evaluating quantities
        # in Fourier space.
        # Note the tolerances are tighter than above
        if guess_mesh_type(psp.r, psp.dr)[1] != "linear"
            rmax = last(psp.r)
            dr_linear = 0.01  # a₀, this is the value used by PseudoDojo

            @testset "Local potential Fourier (on resampled linear mesh) agrees with real" begin
                rmax = cutoff_radius(LocalPotential(), psp)
                r_linear = rmin:dr_linear:rmax
                Vloc_fourier = psp_quantity_evaluator(FourierSpace(), LocalPotential(), psp, r_linear, dr_linear)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(local_potential_integrand(psp, q), rmin, rmax)[1] -
                        4π * psp.Zval / q^2
                    @test ref ≈ Vloc_fourier(q) rtol = 1e-3 atol = 1e-3
                end
            end


            @testset "Nonlocal projector Fouriers (on resampled linear mesh) agree with real" begin
                for l in angular_momenta(psp), n in 1:n_radials(BetaProjector(), psp, l)
                    rmax = cutoff_radius(BetaProjector(), psp, l, n)
                    r_linear = rmin:dr_linear:rmax
                    β_real = psp_quantity_evaluator(RealSpace(), BetaProjector(), psp, l, n)
                    β_fourier = psp_quantity_evaluator(FourierSpace(), BetaProjector(), psp, l, n, r_linear, dr_linear)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, β_real, l, q), rmin, rmax)[1]
                        @test ref ≈ β_fourier(q) rtol = 1e-3 atol = 1e-3
                    end
                end
            end

            @testset "Core charge density Fourier (on resampled linear mesh) agrees with real" begin
                if has_quantity(CoreChargeDensity(), psp)
                    rmax = cutoff_radius(CoreChargeDensity(), psp)
                    r_linear = rmin:dr_linear:rmax
                    ρcore_real = psp_quantity_evaluator(RealSpace(), CoreChargeDensity(), psp)
                    ρcore_fourier = psp_quantity_evaluator(FourierSpace(), CoreChargeDensity(), psp, r_linear, dr_linear)
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, ρcore_real, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρcore_fourier(q) rtol = 1e-5 atol = 1e-5
                    end
                end
            end

            @testset "Pseudo-atomic wavefunction Fouriers (on resampled linear mesh) agree with real" begin
                if has_quantity(ChiProjector(), psp)
                    for l in angular_momenta(psp), n in 1:n_radials(ChiProjector(), psp, l)
                        rmax = cutoff_radius(ChiProjector(), psp, l, n)
                        r_linear = rmin:dr_linear:rmax
                        χ_real = psp_quantity_evaluator(RealSpace(), ChiProjector(), psp, l, n)
                        χ_fourier = psp_quantity_evaluator(FourierSpace(), ChiProjector(), psp, l, n, r_linear, dr_linear)
                        @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                            ref = quadgk(radial_quantity_integrand(psp, χ_real, l, q), rmin, rmax)[1]
                            @test ref ≈ χ_fourier(q) rtol = 1e-3 atol = 1e-3
                        end
                    end
                end
            end

            @testset "Valence charge density Fouriers (on resampled linear mesh) agree with real" begin
                if has_quantity(ValenceChargeDensity(), psp)
                    rmax = cutoff_radius(ValenceChargeDensity(), psp)
                    r_linear = rmin:dr_linear:rmax
                    ρval_real = psp_quantity_evaluator(RealSpace(), ValenceChargeDensity(), psp)
                    ρval_fourier = psp_quantity_evaluator(FourierSpace(), ValenceChargeDensity(), psp, r_linear, dr_linear)
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, ρval_real, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρval_fourier(q) rtol = 1e-5 atol = 1e-5
                    end
                end
            end

            @testset "Pseudo energy correction (on resampled linear mesh)" begin
                q_small = 1e-5
                rmax = cutoff_radius(LocalPotential(), psp)
                r_linear = rmin:dr_linear:rmax
                ref = psp_quantity_evaluator(FourierSpace(), LocalPotential(), psp, r_linear, dr_linear)(q_small) + 4π * psp.Zval / q_small^2
                @test ref ≈ psp_energy_correction(Float64, psp, r_linear, dr_linear) rtol = 5e-5 atol = 5e-5
            end
        end
    end
end
