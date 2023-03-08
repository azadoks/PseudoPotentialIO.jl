@testset "AbstractPsPFile" begin
    for filepath in TEST_FILEPATHS
        file = load_psp_file(filepath)
        @test isa(format(file), AbstractString)
        @test isa(element(file), AbstractString)
        @test -1 <= max_angular_momentum(file) <= 5
        @test 0 <= n_projector_radials(file)
        @test 0 <= n_pseudo_orbital_radials(file)
        @test 0 <= valence_charge(file)  # <= element(file).number
        @test isa(is_norm_conserving(file), Bool)
        @test isa(is_ultrasoft(file), Bool)
        @test isa(is_paw(file), Bool)
        @test formalism(file) in (NormConservingPsP, UltrasoftPsP, ProjectorAugmentedWavePsP)
        @test isa(has_spin_orbit(file), Bool)
        @test relativistic_treatment(file) in (:scalar, :full)
        @test isa(has_core_density(file), Bool)
    end
end
