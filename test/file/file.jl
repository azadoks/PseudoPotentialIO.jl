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
        @test isa(AbstractString, format(file))
        @test isa(AbstractString, element(file))
        @test isa(Bool, has_spin_orbit(file))
        @test relativistic_treatment(file) in (:scalar, :full)
        @test isa(Bool, has_nlcc(file))
        @test 0 <= valence_charge(file)  # <= element(file).number
        @test 0 <= max_angular_momentum(file) <= 5
        @test 0 <= n_projectors(file)
        @test 0 <= n_pseudo_orbitals(file)
        @test isa(Bool, is_norm_conserving(file))
        @test isa(Bool, is_ultrasoft(file))
        @test isa(Bool, is_paw(file))
        @test formalism(file) in (:norm_conserving, :ultrasoft, :paw)
    end
end