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

    for filepath in filepaths
        file = load_psp_file(filepath)
        if (is_norm_conserving(file) | is_ultrasoft(file)) & !has_spin_orbit(file)
            psp = load_psp(file)
            @test isa(AbstractString, element(psp))
            @test isa(Bool, has_spin_orbit(psp))
            @test relativistic_treatment(psp) in (:scalar, :full)
            @test isa(Bool, has_nlcc(psp))
            @test 0 <= valence_charge(file)  # <= element(file).number
            @test 0 <= max_angular_momentum(psp) <= 5
            @test 0 <= n_projectors(psp)
            @test 0 <= n_pseudo_orbitals(psp)
            @test isa(Bool, is_norm_conserving(psp))
            @test isa(Bool, is_ultrasoft(psp))
            @test isa(Bool, is_paw(psp))
            @test formalism(psp) in (:norm_conserving, :ultrasoft, :paw)

            for l in 0:max_angular_momentum(psp)
                @test size(projector_coupling(psp, l)) == (n_projectors(psp, l), n_projectors(psp, l))
                @test eltype(projector_coupling(psp, l)) <: Real
                for n in 1:n_projectors(psp, l)
                    @test projector_coupling(psp, l, n, n) == projector_coupling(psp, l, n)
                    for m in (n + 1):n_projectors(psp, l)
                        @test isa(Real, projector_coupling(psp, l, n, m))
                        @test projector_coupling(psp, l, n, m) == projector_coupling(psp, l, m, n)
                    end
                end
            end
        end
    end
end
