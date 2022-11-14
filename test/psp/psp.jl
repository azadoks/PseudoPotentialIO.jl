#TODO hide deprecation message coming from Quaternions.jl through Rotations.jl
@testset "AbstractPsP" begin
    filepaths = []
    for dir in families
        (_, _, files) = first(walkdir(dir))
        for file in files
            if file[1] != '.'  # Hack to avoid hidden files
                push!(filepaths, joinpath(dir, file))
            end
        end
    end

    @testset "$(splitpath(filepath)[end])" for filepath in filepaths
        file = load_psp_file(filepath)
        if (is_norm_conserving(file) | is_ultrasoft(file)) & !has_spin_orbit(file)
            psp = load_psp(file)
            @test isa(element(psp), AbstractString)
            @test -1 <= max_angular_momentum(psp) <= 5
            @test 0 <= n_projectors(psp)
            @test 0 <= n_pseudo_orbitals(psp)
            @test 0 <= valence_charge(file)  # <= element(file).number
            @test isa(is_norm_conserving(psp), Bool)
            @test isa(is_ultrasoft(psp), Bool)
            @test isa(is_paw(psp), Bool)
            @test isa(has_spin_orbit(psp), Bool)
            @test isa(has_nlcc(psp), Bool)

            for l in 0:max_angular_momentum(psp)
                @test size(projector_coupling(psp, l)) ==
                      (n_projectors(psp, l), n_projectors(psp, l))
                @test eltype(projector_coupling(psp, l)) <: Real
                for n in 1:n_projectors(psp, l)
                    @test projector_coupling(psp, l, n, n) == projector_coupling(psp, l, n)
                    for m in (n + 1):n_projectors(psp, l)
                        @test isa(projector_coupling(psp, l, n, m), Real)
                        @test projector_coupling(psp, l, n, m) ==
                              projector_coupling(psp, l, m, n)
                    end
                end
            end

            R = rand(3)
            R_rot = QuatRotation(rand(RotMatrix{3})) * R

            @test local_potential_real(psp, R) ≈ local_potential_real(psp, R_rot)
            for l in 0:max_angular_momentum(psp), n in 1:n_projectors(psp, l)
                @test projector_real(psp, l, n, R) ≈ projector_real(psp, l, n, R_rot)
            end
            for l in 0:max_angular_momentum(psp), n in 1:n_pseudo_orbitals(psp, l)
                @test pseudo_orbital_real(psp, l, n, R) ≈
                      pseudo_orbital_real(psp, l, n, R_rot)
            end
            ρval_R = valence_charge_density_real(psp, R)
            @test typeof(ρval_R) <: Union{Nothing,Real}
            if !isnothing(ρval_R)
                @test ρval_R ≈ valence_charge_density_real(psp, R_rot)
            end
            ρcore_R = core_charge_density_real(psp, R)
            @test typeof(ρcore_R) <: Union{Nothing,Real}
            if !isnothing(ρcore_R)
                @test ρcore_R ≈ core_charge_density_real(psp, R_rot)
            end

            K = rand(3) .+ eps(Float64)
            K_rot = QuatRotation(rand(RotMatrix{3})) * K

            @test local_potential_fourier(psp, K) ≈ local_potential_fourier(psp, K_rot)
            for l in 0:max_angular_momentum(psp), n in 1:n_projectors(psp, l)
                @test projector_fourier(psp, l, n, K) ≈ projector_fourier(psp, l, n, K_rot)
            end
            for l in 0:max_angular_momentum(psp), n in 1:n_pseudo_orbitals(psp, l)
                @test pseudo_orbital_fourier(psp, l, n, K) ≈
                      pseudo_orbital_fourier(psp, l, n, K_rot)
            end
            ρval_K = valence_charge_density_fourier(psp, K)
            @test typeof(ρval_K) <: Union{Nothing,Real}
            if !isnothing(ρval_K)
                @test ρval_K ≈ valence_charge_density_fourier(psp, K_rot)
            end
            ρcore_K = core_charge_density_fourier(psp, K)
            @test typeof(ρcore_K) <: Union{Nothing,Real}
            if !isnothing(ρcore_K)
                @test ρcore_K ≈ core_charge_density_fourier(psp, K_rot)
            end

            @test relativistic_treatment(psp) in (:scalar, :full)
            @test formalism(psp) in (:norm_conserving, :ultrasoft, :paw)
        end
    end
end
