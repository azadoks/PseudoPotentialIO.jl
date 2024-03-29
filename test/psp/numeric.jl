import PseudoPotentialIO: build_interpolator_real

#TODO work on tightening tolerances
@testset "Numeric" begin
    function local_potential_integrand(psp::NumericPsP, q)
        Vloc = local_potential_real(psp)
        return r -> 4π * r * sphericalbesselj(0, q * r) * (r * Vloc(r) + psp.Zval)
    end

    function radial_quantity_integrand(::NumericPsP, f, l::Int, q)
        return r -> 4π * sphericalbesselj(l, q * r) * f(r)
    end

    @testset "$(name)" for (name, filepath) in NUMERIC_CASE_FILEPATHS
        @testset "Truncation tol=$(tol)" for tol in [nothing, 1e-10, 1e-8, 1e-6]
            psp = load_psp(filepath)
            rmin = max(first(psp.r), eps(eltype(psp.r)))

            @testset "Local potential Fourier agrees with real" begin
                rmax = local_potential_cutoff_radius(psp)
                Vloc_fourier = local_potential_fourier(psp; tol)
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(local_potential_integrand(psp, q), rmin, rmax)[1] -
                        4π * psp.Zval / q^2
                    @test ref ≈ Vloc_fourier(q) rtol = 1e-2 atol = 1e-2
                end
            end

            @testset "Nonlocal projector Fouriers agree with real" begin
                for l in 0:(psp.lmax), n in eachindex(psp.β[l])
                    itp = build_interpolator_real(psp.β[l][n], psp.r)
                    rmax = projector_cutoff_radius(psp, l, n)
                    β_fourier = projector_fourier(psp, l, n; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(radial_quantity_integrand(psp, itp, l, q), rmin, rmax)[1]
                        @test ref ≈ β_fourier(q) rtol = 1e-1 atol = 1e-1
                    end
                end
            end

            @testset "Core charge density Fourier agrees with real" begin
                if !isnothing(psp.ρcore)
                    itp = build_interpolator_real(psp.ρcore, psp.r)
                    rmax = core_charge_density_cutoff_radius(psp)
                    ρcore_fourier = core_charge_density_fourier(psp; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(radial_quantity_integrand(psp, itp, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρcore_fourier(q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end

            @testset "Pseudo energy correction" begin
                q_small = 1e-5
                ref = local_potential_fourier(psp)(q_small) + 4π * psp.Zval / q_small^2
                @test ref ≈ pseudo_energy_correction(Float64, psp) rtol = 1e-2 atol = 1e-2
            end

            @testset "Pseudo-atomic wavefunction Fouriers agree with real" begin
                if !isnothing(psp.χ)
                    for l in 0:(psp.lmax), n in eachindex(psp.χ[l])
                        itp = build_interpolator_real(psp.χ[l][n], psp.r)
                        rmax = chi_function_cutoff_radius(psp, l, n)
                        χ_fourier = chi_function_fourier(psp, l, n; tol)
                        for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                            ref = quadgk(radial_quantity_integrand(psp, itp, l, q), rmin, rmax)[1]
                            @test ref ≈ χ_fourier(q) rtol = 1e-1 atol = 1e-1
                        end
                    end
                end
            end

            @testset "Valence charge density Fouriers agree with real" begin
                if !isnothing(psp.ρval)
                    itp = build_interpolator_real(psp.ρval, psp.r)
                    rmax = valence_charge_density_cutoff_radius(psp)
                    ρval_fourier = valence_charge_density_fourier(psp; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(radial_quantity_integrand(psp, itp, 0, q), rmin, rmax)[1]
                        @test ref ≈ ρval_fourier(q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end
        end
    end
end
