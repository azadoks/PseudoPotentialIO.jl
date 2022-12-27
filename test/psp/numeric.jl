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

    function local_potential_integral(psp::NumericPsP, itp, q::T, r::T)::T where {T<:Real}
        return 4π * r * fast_sphericalbesselj0(q * r) * (r * itp(r) + psp.Zval)
    end

    function wavefunction_like_integral(::NumericPsP, itp, l::Int, q::T,
                                        r::T)::T where {T<:Real}
        return 4π * r^2 * fast_sphericalbesselj(l, q * r) * itp(r)
    end

    function charge_density_integral(::NumericPsP, itp, q::T, r::T)::T where {T<:Real}
        return 4π * r^2 * fast_sphericalbesselj0(q * r) * itp(r)
    end

    @testset "$(splitpath(filepath)[end])" for filepath in filepaths
        psp = load_psp(filepath)
        rmin = first(psp.r)

        @testset "Local potential Fourier agrees with real" begin
            itp = build_interpolator(psp.Vloc, psp.r)
            rmax = psp.r[lastindex(psp.Vloc)]
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                ref = quadgk(r -> local_potential_integral(psp, itp, q, r), rmin, rmax)[1] -
                      4π * psp.Zval / q^2
                @test ref ≈ local_potential_fourier(psp, q) rtol = 1e-2 atol = 1e-2
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in 0:(psp.lmax), n in eachindex(psp.β[l])
                itp = build_interpolator(psp.β[l][n], psp.r)
                rmax = psp.r[lastindex(psp.β[l][n])]
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> wavefunction_like_integral(psp, itp, l, q, r), rmin,
                                 rmax)[1]
                    @test ref ≈ projector_fourier(psp, l, n, q) rtol = 1e-2 atol = 1e-2
                end
            end
        end

        @testset "Core charge density Fourier agrees with real" begin
            if !isnothing(psp.ρcore)
                itp = build_interpolator(psp.ρcore, psp.r)
                rmax = psp.r[lastindex(psp.ρcore)]
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> charge_density_integral(psp, itp, q, r), rmin, rmax)[1]
                    @test ref ≈ core_charge_density_fourier(psp, q) rtol = 1e-2 atol = 1e-2
                end
            end
        end

        @testset "Pseudo energy correction" begin
            q_small = 1e-5
            ref = local_potential_fourier(psp, q_small) + 4π * psp.Zval / q_small^2
            @test ref ≈ pseudo_energy_correction(psp) rtol = 1e-2 atol = 1e-2
        end

        @testset "Pseudo-atomic wavefunction Fouriers agree with real" begin
            if !isnothing(psp.ϕ̃)
                for l in 0:(psp.lmax), n in eachindex(psp.ϕ̃[l])
                    itp = build_interpolator(psp.ϕ̃[l][n], psp.r)
                    rmax = psp.r[lastindex(psp.ϕ̃[l][n])]
                    for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                        ref = quadgk(r -> wavefunction_like_integral(psp, itp, l, q, r),
                                     rmin, rmax)[1]
                        @test ref ≈ pseudo_orbital_fourier(psp, l, n, q) rtol = 1e-2 atol = 1e-2
                    end
                end
            end
        end

        @testset "Valence charge density Fouriers agree with real" begin
            if !isnothing(psp.ρval)
                itp = build_interpolator(psp.ρval, psp.r)
                rmax = psp.r[lastindex(psp.ρval)]
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> charge_density_integral(psp, itp, q, r), rmin, rmax)[1]
                    @test ref ≈ valence_charge_density_fourier(psp, q) rtol = 1e-2 atol = 1e-2
                end
            end
        end
    end
end
