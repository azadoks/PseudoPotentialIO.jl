import PseudoPotentialIO: linear_mesh, logarithmic_mesh1, logarithmic_mesh2,
                          guess_mesh_type, sphericalbesselj_fast, trapezoid, simpson
@testset "common" begin
    @testset "guess_mesh_type" begin
        @testset "[linear] Mg_nc-fr-04_pbesol_stringent.upf" begin
            psp = UpfPsP("upf2/Mg_nc-fr-04_pbesol_stringent.upf")
            type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
            @test type_guess == "linear"
            @test a_guess ≈ 0.01
            @test b_guess ≈ -0.01
        end

        @testset "[logarithmic type 1] Al.pbe-n-kjpaw_psl.1.0.0.upf" begin
            psp = UpfPsP("upf2/Al.pbe-n-kjpaw_psl.1.0.0.upf")
            type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
            @test type_guess == "log_1"
            @test a_guess ≈ psp.mesh.dx
            @test b_guess ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
        end

        @testset "[logarithmic type 1] he-q2.upf" begin
            psp = UpfPsP("upf2/he-q2.upf")
            type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
            @test type_guess == "log_1"
            @test a_guess ≈ psp.mesh.dx
            @test b_guess ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
        end

        @testset "[logarithmic type 2] constructed" begin
            #TODO: find a log_2 type mesh example
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

    @testset "sphericalbesselj_fast" begin
        for l in 0:5
            for x in (rand(100) .* 100)
                @test sphericalbesselj(l, x) ≈ sphericalbesselj_fast(l, x) atol = 1e-8
            end
        end
        for l in 6:10
            @test_throws "The case l = $l is not implemented" sphericalbesselj_fast(l, 1.0)
            @test sphericalbesselj_fast(l, 0.0) == 0.0
        end
    end

    # f(x) = exp(-x) * sin(x)
    # @testset "trapezoid" begin
    #     @testset "linear mesh" begin
    #         i = 1:101
    #         x = linear_mesh.(i, 0.01, -0.01)
    #         dx = 0.01
    #         fx = f.(x)
    #         reference = quadgk(cubic_spline_interpolation(minimum(x):dx:maximum(x), fx),
    #                            first(x), last(x))[1]
    #         @test trapezoid(fx, dx) ≈ reference rtol = 1e-4
    #     end

    #     @testset "logarithmic mesh" begin
    #         i = 1:101
    #         a = 0.01
    #         b = exp(-7.0) / 0.66
    #         x = logarithmic_mesh1.(i, a, b)
    #         dx = a .* b .* exp.((i .- 1) .* a)
    #         fx = f.(x)
    #         reference = quadgk(linear_interpolation((x,), fx), first(x), last(x))[1]
    #         @test trapezoid(fx, dx) ≈ reference rtol = 1e-4
    #     end
    # end

    # @testset "simpson" begin
    #     @testset "odd linear mesh" begin
    #         i = 1:101
    #         x = linear_mesh.(i, 0.01, -0.01)
    #         dx = 0.01
    #         fx = f.(x)
    #         reference = quadgk(cubic_spline_interpolation(minimum(x):dx:maximum(x), fx),
    #                            first(x), last(x))[1]
    #         @test simpson(fx, dx) ≈ reference rtol = 1e-4
    #     end

    #     @testset "odd logarithmic mesh" begin
    #         i = 1:101
    #         a = 0.01
    #         b = exp(-7.0) / 0.66
    #         x = logarithmic_mesh1.(i, a, b)
    #         dx = a .* b .* exp.((i .- 1) .* a)
    #         fx = f.(x)
    #         reference = quadgk(linear_interpolation((x,), fx), first(x), last(x))[1]
    #         @test simpson(fx, dx) ≈ reference rtol = 1e-4
    #     end

    #     @testset "even linear mesh" begin
    #         i = 1:100
    #         x = linear_mesh.(i, 0.01, -0.01)
    #         dx = 0.01
    #         fx = f.(x)
    #         reference = quadgk(cubic_spline_interpolation(minimum(x):dx:maximum(x), fx),
    #                            first(x), last(x))[1]
    #         @test simpson(fx, dx) ≈ reference rtol = 1e-4
    #     end

    #     @testset "even logarithmic mesh" begin
    #         i = 1:100
    #         a = 0.01
    #         b = exp(-7.0) / 0.66
    #         x = logarithmic_mesh1.(i, a, b)
    #         dx = a .* b .* exp.((i .- 1) .* a)
    #         fx = f.(x)
    #         reference = quadgk(linear_interpolation((x,), fx), first(x), last(x))[1]
    #         @test simpson(fx, dx) ≈ reference rtol = 1e-4
    #     end
    # end
end