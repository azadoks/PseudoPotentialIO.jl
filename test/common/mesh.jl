import PseudoPotentialIO: linear_mesh, logarithmic_mesh1, logarithmic_mesh2, guess_mesh_type

@testset "Mesh type guess" begin
    @testset "[linear] Mg.upf" begin
        psp = load_psp_file(UPF2_CASE_FILEPATHS["Mg.upf"])
        type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
        @test type_guess == "linear"
        @test a_guess ≈ 0.01
        @test b_guess ≈ -0.01
    end

    @testset "[logarithmic type 1] Al.pbe-n-kjpaw_psl.1.0.0.UPF" begin
        psp = load_psp_file(UPF2_CASE_FILEPATHS["Al.pbe-n-kjpaw_psl.1.0.0.UPF"])
        type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
        @test type_guess == "log_1"
        @test a_guess ≈ psp.mesh.dx
        @test b_guess ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
    end

    @testset "[logarithmic type 1] He.pbe-hgh.UPF" begin
        psp = load_psp_file(UPF2_CASE_FILEPATHS["He.pbe-hgh.UPF"])
        type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
        @test type_guess == "log_1"
        @test a_guess ≈ psp.mesh.dx
        @test b_guess ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
    end

    @testset "[logarithmic type 2] constructed" begin
        #TODO: find a log_2 type mesh example file
        i = 1:1000
        a = 0.0125  # dx
        b = exp(-7.0) / 0.66  # exp(xmin) / zmesh
        log2_mesh = logarithmic_mesh2.(i, a, b)
        log2_dmesh = a .* b .* exp.((i .- 1) .* a)  # Derivative of log2 mesh w.r.t. `i`

        type_guess, a_guess, b_guess = guess_mesh_type(log2_mesh, log2_dmesh)
        @test type_guess == "log_2"
        @test a_guess ≈ a
        @test b_guess ≈ b
    end
end
