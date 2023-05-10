import PseudoPotentialIO: guess_mesh_type

#TODO work on tightening tolerances

#* Note: q-values are equivalent to the largest q at an Ecut via q^2 = 2Ecut
#*       Ecut up to ≈250 Ha is considered.
@testset "Numeric" begin
    function local_potential_integrand(psp::NumericPsP, q)
        Vloc = psp_quantity_evaluator(psp, LocalPotential(), RealSpace())
        return r -> 4π * r * sphericalbesselj(0, q * r) * (r * Vloc(r) + psp.Zval)
    end

    function radial_quantity_integrand(::NumericPsP, f, l::Int, q)
        return r -> 4π * sphericalbesselj(l, q * r) * f(r)
    end

    @testset "$(name)" for (name, filepath) in NUMERIC_CASE_FILEPATHS
        psp = load_psp(filepath)
        rmin = max(first(psp.r), eps(eltype(psp.r)))

        @testset "Local potential Fourier agrees with real" begin
            rmax = cutoff_radius(psp, LocalPotential())
            Vloc_fourier = psp_quantity_evaluator(psp, LocalPotential(), FourierSpace())
            @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                ref = quadgk(local_potential_integrand(psp, q), rmin, rmax)[1] -
                    4π * psp.Zval / q^2
                @test ref ≈ Vloc_fourier(q) rtol = 1e-2 atol = 1e-2
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in angular_momenta(psp), n in 1:n_radials(psp, BetaProjector(), l)
                rmax = psp.r[lastindex(get_quantity(psp, BetaProjector(), l, n))]
                β_real = psp_quantity_evaluator(psp, BetaProjector(), l, n, RealSpace())
                β_fourier = psp_quantity_evaluator(psp, BetaProjector(), l, n, FourierSpace())
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, β_real, l, q), rmin, rmax)[1]
                    @test ref ≈ β_fourier(q) rtol = 1e-2 atol = 1e-2
                end
            end
        end

        @testset "Core charge density Fourier agrees with real" begin
            if has_quantity(psp, CoreChargeDensity())
                rmax = cutoff_radius(psp, CoreChargeDensity())
                ρcore_real = psp_quantity_evaluator(psp, CoreChargeDensity(), RealSpace())
                ρcore_fourier = psp_quantity_evaluator(psp, CoreChargeDensity(), FourierSpace())
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, ρcore_real, 0, q), rmin, rmax)[1]
                    @test ref ≈ ρcore_fourier(q) rtol = 1e-4 atol = 1e-4
                end
            end
        end

        @testset "Pseudo-atomic wavefunction Fouriers agree with real" begin
            if has_quantity(psp, ChiProjector())
                for l in angular_momenta(psp), n in 1:n_radials(psp, ChiProjector(), l)
                    rmax = psp.r[lastindex(get_quantity(psp, ChiProjector(), l, n))]
                    χ_real = psp_quantity_evaluator(psp, ChiProjector(), l, n, RealSpace())
                    χ_fourier = psp_quantity_evaluator(psp, ChiProjector(), l, n, FourierSpace())
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, χ_real, l, q), rmin, rmax)[1]
                        @test ref ≈ χ_fourier(q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end
        end

        @testset "Valence charge density Fouriers agree with real" begin
            if has_quantity(psp, ValenceChargeDensity())
                rmax = cutoff_radius(psp, ValenceChargeDensity())
                ρval_real = psp_quantity_evaluator(psp, ValenceChargeDensity(), RealSpace())
                ρval_fourier = psp_quantity_evaluator(psp, ValenceChargeDensity(), FourierSpace())
                @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(radial_quantity_integrand(psp, ρval_real, 0, q), rmin, rmax)[1]
                    @test ref ≈ ρval_fourier(q) rtol = 1e-4 atol = 1e-4
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = psp_quantity_evaluator(psp, LocalPotential(), FourierSpace())(q_small) + 4π * psp.Zval / q_small^2
            @test ref ≈ psp_energy_correction(Float64, psp) rtol = 1e-4 atol = 1e-4
        end

        # Test resampling log-mesh pseudopotentials on a linear mesh using cubic splines
        # before evaluating quantities in Fourier space. Empirically, this improves
        # the quality of the Fourier quantities and at high densities will cause
        # the quality of the Trapezoid integrator to approach the Simpson integrator.
        # Note the tolerances are tighter than above.
        if guess_mesh_type(psp.r, psp.dr)[1] != "linear"
            rmax = last(psp.r)
            Δr_linear = 0.001  # a₀, this is 1/10 the value used by PseudoDojo

            @testset "Local potential Fourier (on resampled linear mesh) agrees with real" begin
                rmax = cutoff_radius(psp, LocalPotential())
                r_linear = rmin:Δr_linear:rmax
                dr_linear = fill(Δr_linear, length(r_linear))
                Vloc_fourier = psp_quantity_evaluator(psp, LocalPotential(), FourierSpace(), r_linear, dr_linear, Δr_linear)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                    ref = quadgk(local_potential_integrand(psp, q), rmin, rmax)[1] -
                        4π * psp.Zval / q^2
                    @test ref ≈ Vloc_fourier(q) rtol = 1e-3 atol = 1e-3
                end
            end


            @testset "Nonlocal projector Fouriers (on resampled linear mesh) agree with real" begin
                for l in angular_momenta(psp), n in 1:n_radials(psp, BetaProjector(), l)
                    rmax = psp.r[lastindex(get_quantity(psp, BetaProjector(), l, n))]
                    r_linear = rmin:Δr_linear:rmax
                    dr_linear = fill(Δr_linear, length(r_linear))
                    β_real = psp_quantity_evaluator(psp, BetaProjector(), l, n, RealSpace())
                    β_fourier = psp_quantity_evaluator(psp, BetaProjector(), l, n, FourierSpace(), r_linear, dr_linear, Δr_linear)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, β_real, l, q), rmin, rmax)[1]
                        @test ref ≈ β_fourier(q) rtol = 1e-3 atol = 1e-3
                    end
                end
            end

            @testset "Core charge density Fourier (on resampled linear mesh) agrees with real" begin
                if has_quantity(psp, CoreChargeDensity())
                    rmax = cutoff_radius(psp, CoreChargeDensity())
                    r_linear = rmin:Δr_linear:rmax
                    dr_linear = fill(Δr_linear, length(r_linear))
                    ρcore_real = psp_quantity_evaluator(psp, CoreChargeDensity(), RealSpace())
                    ρcore_fourier = psp_quantity_evaluator(psp, CoreChargeDensity(), FourierSpace(), r_linear, dr_linear, Δr_linear)
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, ρcore_real, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρcore_fourier(q) rtol = 1e-5 atol = 1e-5
                    end
                end
            end

            @testset "Pseudo-atomic wavefunction Fouriers (on resampled linear mesh) agree with real" begin
                if has_quantity(psp, ChiProjector())
                    for l in angular_momenta(psp), n in 1:n_radials(psp, ChiProjector(), l)
                        rmax = psp.r[lastindex(get_quantity(psp, ChiProjector(), l, n))]
                        r_linear = rmin:Δr_linear:rmax
                        dr_linear = fill(Δr_linear, length(r_linear))
                        χ_real = psp_quantity_evaluator(psp, ChiProjector(), l, n, RealSpace())
                        χ_fourier = psp_quantity_evaluator(psp, ChiProjector(), l, n, FourierSpace(), r_linear, dr_linear, Δr_linear)
                        @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                            ref = quadgk(radial_quantity_integrand(psp, χ_real, l, q), rmin, rmax)[1]
                            @test ref ≈ χ_fourier(q) rtol = 1e-3 atol = 1e-3
                        end
                    end
                end
            end

            @testset "Valence charge density Fouriers (on resampled linear mesh) agree with real" begin
                if has_quantity(psp, ValenceChargeDensity())
                    rmax = cutoff_radius(psp, ValenceChargeDensity())
                    r_linear = rmin:Δr_linear:rmax
                    dr_linear = fill(Δr_linear, length(r_linear))
                    ρval_real = psp_quantity_evaluator(psp, ValenceChargeDensity(), RealSpace())
                    ρval_fourier = psp_quantity_evaluator(psp, ValenceChargeDensity(), FourierSpace(), r_linear, dr_linear, Δr_linear)
                    @testset "$(q)" for q in (0.01, 0.5, 2.5, 5.0, 10.0, 22.0)
                        ref = quadgk(radial_quantity_integrand(psp, ρval_real, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρval_fourier(q) rtol = 1e-5 atol = 1e-5
                    end
                end
            end

            @testset "Pseudo energy correction (on resampled linear mesh)" begin
                q_small = 1e-5
                rmax = cutoff_radius(psp, LocalPotential())
                r_linear = rmin:Δr_linear:rmax
                dr_linear = fill(Δr_linear, length(r_linear))
                ref = psp_quantity_evaluator(psp, LocalPotential(), FourierSpace(), r_linear, dr_linear, Δr_linear)(q_small) + 4π * psp.Zval / q_small^2
                @test ref ≈ psp_energy_correction(Float64, psp, r_linear, dr_linear, Δr_linear) rtol = 5e-5 atol = 5e-5
            end
        end
    end
end
