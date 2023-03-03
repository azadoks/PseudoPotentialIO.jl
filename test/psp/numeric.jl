import PseudoPotentialIO: fast_sphericalbesselj, fast_sphericalbesselj0
import PseudoPotentialIO: build_interpolator_real

#TODO work on tightening tolerances
@testset "Numeric" begin
    filepaths = [upf1_filepaths["B_pbe_v1.01.uspp.F.UPF"],
                 upf1_filepaths["si_pbesol_v1.uspp.F.UPF"],
                 upf2_filepaths["He.pbe-hgh.UPF"],
                 upf2_filepaths["Al.pbe-hgh.UPF"],
                 upf2_filepaths["Si.pbe-n-rrkjus_psl.1.0.0.UPF"],
                 upf2_filepaths["H.upf"],
                 upf2_filepaths["Fe.upf"],
                 psp8_filepaths["H.psp8"],
                 psp8_filepaths["Fe.psp8"]]

    function local_potential_integrand(psp::NumericPsP, q)
        Vloc = local_potential_real(psp)
        return r -> 4π * r * fast_sphericalbesselj0(q * r) * (r * Vloc(r) + psp.Zval)
    end

    function wavefunction_like_integrand(::NumericPsP, itp, l::Int, q)
        jₗ = fast_sphericalbesselj(l)
        return r -> 4π * jₗ(q * r) * itp(r)
    end

    function charge_density_integrand(::NumericPsP, itp, q)
        return r -> 4π * fast_sphericalbesselj0(q * r) * itp(r)
    end

    @testset "$(splitpath(filepath)[end])" for filepath in filepaths
        @testset "Truncation tol=$(tol)" for tol in [nothing, 1e-10, 1e-8, 1e-6]
            psp = load_psp(filepath)
            rmin = max(first(psp.r), eps(eltype(psp.r)))

            @testset "Local potential Fourier agrees with real" begin
                rmax = psp.r[lastindex(psp.Vloc)]
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
                    rmax = psp.r[lastindex(psp.β[l][n])]
                    β_fourier = projector_fourier(psp, l, n; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(wavefunction_like_integrand(psp, itp, l, q), rmin, rmax)[1]
                        @test ref ≈ β_fourier(q) rtol = 1e-1 atol = 1e-1
                    end
                end
            end

            @testset "Core charge density Fourier agrees with real" begin
                if !isnothing(psp.ρcore)
                    itp = build_interpolator_real(psp.ρcore, psp.r)
                    rmax = psp.r[lastindex(psp.ρcore)]
                    ρcore_fourier = core_charge_density_fourier(psp; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(charge_density_integrand(psp, itp, q), rmin, rmax)[1]
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
                if !isnothing(psp.ϕ)
                    for l in 0:(psp.lmax), n in eachindex(psp.ϕ[l])
                        itp = build_interpolator_real(psp.ϕ[l][n], psp.r)
                        rmax = psp.r[lastindex(psp.ϕ[l][n])]
                        ϕ_fourier = pseudo_orbital_fourier(psp, l, n; tol)
                        for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                            ref = quadgk(wavefunction_like_integrand(psp, itp, l, q), rmin, rmax)[1]
                            @test ref ≈ ϕ_fourier(q) rtol = 1e-1 atol = 1e-1
                        end
                    end
                end
            end

            @testset "Valence charge density Fouriers agree with real" begin
                if !isnothing(psp.ρval)
                    itp = build_interpolator_real(psp.ρval, psp.r)
                    rmax = psp.r[lastindex(psp.ρval)]
                    ρval_fourier = valence_charge_density_fourier(psp; tol)
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(charge_density_integrand(psp, itp, q), rmin, rmax)[1]
                        @test ref ≈ ρval_fourier(q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end
        end
    end
end
