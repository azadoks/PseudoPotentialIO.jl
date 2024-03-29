function random_versor()
    r, i, j, k = randn(4)
    quat_inorm = inv(sqrt(r * r + i * i + j * j + k * k))
    return quat(r * quat_inorm, i * quat_inorm, j * quat_inorm, k * quat_inorm)
end

function rotate_vector(q::Quaternion, u::AbstractVector)
    q_u = Quaternion(0, u[1], u[2], u[3])
    q_v = q * q_u * conj(q)
    return [imag_part(q_v)...]
end

@testset "AbstractPsP" begin
    @testset "$(splitpath(filepath)[end])" for filepath in TEST_FILEPATHS
        file = load_psp_file(filepath)
        if (is_norm_conserving(file) | is_ultrasoft(file)) & !has_spin_orbit(file)
            psp = load_psp(file)
            @test isa(element(psp), AbstractString)
            @test -1 <= max_angular_momentum(psp) <= 5
            @test 0 <= n_projector_radials(psp)
            @test 0 <= n_chi_function_radials(psp)
            @test 0 <= valence_charge(file)  # <= element(file).number
            @test isa(is_norm_conserving(psp), Bool)
            @test isa(is_ultrasoft(psp), Bool)
            @test isa(is_paw(psp), Bool)
            @test isa(has_spin_orbit(psp), Bool)
            @test isa(has_core_density(psp), Bool)

            for l in angular_momenta(psp)
                @test size(projector_coupling(psp, l)) ==
                      (n_projector_radials(psp, l), n_projector_radials(psp, l))
                @test eltype(projector_coupling(psp, l)) <: Real
                for n in projector_radial_indices(psp, l)
                    @test projector_coupling(psp, l, n, n) == projector_coupling(psp, l, n)
                    for m in (n + 1):n_projector_radials(psp, l)
                        @test isa(projector_coupling(psp, l, n, m), Real)
                        @test projector_coupling(psp, l, n, m) ==
                              projector_coupling(psp, l, m, n)
                    end
                end
            end

            R = rand(3)
            r = norm(R)
            R̂ = R ./ norm(r)
            if isfinite(pseudo_cutoff_radius(psp))
                r = r / sqrt(3) * pseudo_cutoff_radius(psp)
            end
            R = R̂ * r
            R_rot = rotate_vector(random_versor(), R)

            @test local_potential_real(psp)(R) ≈ local_potential_real(psp)(R_rot)
            for l in angular_momenta(psp), n in projector_radial_indices(psp, l)
                @test projector_real(psp, l, n)(R) ≈ projector_real(psp, l, n)(R_rot)
            end
            for l in angular_momenta(psp), n in chi_function_radial_indices(psp, l)
                @test chi_function_real(psp, l, n)(R) ≈
                      chi_function_real(psp, l, n)(R_rot)
            end
            ρval_R = valence_charge_density_real(psp)(R)
            @test typeof(ρval_R) <: Union{Nothing,Real}
            if !isnothing(ρval_R)
                @test ρval_R ≈ valence_charge_density_real(psp)(R_rot)
            end
            ρcore_R = core_charge_density_real(psp)(R)
            @test typeof(ρcore_R) <: Union{Nothing,Real}
            if !isnothing(ρcore_R)
                @test ρcore_R ≈ core_charge_density_real(psp)(R_rot)
            end

            K = rand(3) .+ eps(Float64)
            K_rot = rotate_vector(random_versor(), K)

            @test local_potential_fourier(psp)(K) ≈ local_potential_fourier(psp)(K_rot)
            for l in angular_momenta(psp), n in projector_radial_indices(psp, l)
                @test projector_fourier(psp, l, n)(K) ≈ projector_fourier(psp, l, n)(K_rot)
            end
            for l in angular_momenta(psp), n in chi_function_radial_indices(psp, l)
                @test chi_function_fourier(psp, l, n)(K) ≈
                      chi_function_fourier(psp, l, n)(K_rot)
            end
            ρval_K = valence_charge_density_fourier(psp)(K)
            @test typeof(ρval_K) <: Union{Nothing,Real}
            if !isnothing(ρval_K)
                @test ρval_K ≈ valence_charge_density_fourier(psp)(K_rot)
            end
            ρcore_K = core_charge_density_fourier(psp)(K)
            @test typeof(ρcore_K) <: Union{Nothing,Real}
            if !isnothing(ρcore_K)
                @test ρcore_K ≈ core_charge_density_fourier(psp)(K_rot)
            end

            @test relativistic_treatment(psp) in (:scalar, :full)
            @test formalism(psp) in
                  (NormConservingPsP, UltrasoftPsP, ProjectorAugmentedWavePsP)
        end
    end
end
