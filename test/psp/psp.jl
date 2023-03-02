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

            @test relativistic_treatment(psp) in (:scalar, :full)
            @test formalism(psp) in (NormConservingPsP, UltrasoftPsP, ProjectorAugmentedWavePsP)
        end
    end
end
