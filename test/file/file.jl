@testset "AbstractPsPFile" begin
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
        @test isa(format(file), AbstractString)
        @test isa(element(file), AbstractString)
        @test -1 <= max_angular_momentum(file) <= 5
        @test 0 <= n_projectors(file)
        @test 0 <= n_pseudo_orbitals(file)
        @test 0 <= valence_charge(file)  # <= element(file).number
        @test isa(is_norm_conserving(file), Bool)
        @test isa(is_ultrasoft(file), Bool)
        @test isa(is_paw(file), Bool)
        @test formalism(file) in (:norm_conserving, :ultrasoft, :paw)
        @test isa(has_spin_orbit(file), Bool)
        @test relativistic_treatment(file) in (:scalar, :full)
        @test isa(has_nlcc(file), Bool)
    end
end