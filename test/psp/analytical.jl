@testset "Analytical" begin
    dirs = [artifact"hgh_lda_hgh", artifact"hgh_pbe_hgh"]

    filepaths = []
    for dir in dirs
        (_, _, files) = first(walkdir(dir))
        for file in files
            if file[1] != '.'  # Hack to avoid hidden files
                push!(filepaths, joinpath(dir, file))
                #TODO remove
                break
                #TODO
            end
        end
    end

    function local_potential_integral(psp, q, r)
        return 4π * r * sphericalbesselj(0, q * r) *
               (r * local_potential_real(psp, r) + valence_charge(psp))
    end

    function projector_integral(psp, l, n, q, r)
        return 4π * r^2 * sphericalbesselj(l, q * r) * projector_real(psp, l, n, r)
    end

    @testset "$(splitpath(filepath)[end])" for filepath in filepaths
        psp = load_psp(filepath)

        @testset "Local potential Fourier agrees with real" begin
            for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                ref = quadgk(r -> local_potential_integral(psp, q, r), 0, Inf)[1] -
                      4π * psp.Zval / q^2
                @test ref ≈ local_potential_fourier(psp, q) rtol = 1e-3 atol = 1e-3
            end
        end

        @testset "Nonlocal projector Fouriers agree with real" begin
            for l in 0:(psp.lmax), n in eachindex(psp.β[l])
                for q in (0.01, 0.5, 2.5, 5.0, 10.0, 50.0)
                    ref = quadgk(r -> projector_integral(psp, l, n, q, r), 0, Inf)[1]
                    @test ref ≈ projector_fourier(psp, l, n, q) rtol = 1e-3 atol = 1e-3
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
