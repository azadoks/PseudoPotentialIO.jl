import PseudoPotentialIO: linear_mesh, logarithmic_mesh1, logarithmic_mesh2
import PseudoPotentialIO: guess_mesh_type, fast_sphericalbesselj, fast_sphericalbesselj0
import PseudoPotentialIO: trapezoid, simpson

@testset "common" begin
    @testset "guess_mesh_type" begin
        @testset "[linear] Mg.upf" begin
            psp = load_psp_file(upf2_filepaths["Mg.upf"])
            type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
            @test type_guess == "linear"
            @test a_guess ≈ 0.01
            @test b_guess ≈ -0.01
        end

        @testset "[logarithmic type 1] Al.pbe-n-kjpaw_psl.1.0.0.UPF" begin
            psp = load_psp_file(upf2_filepaths["Al.pbe-n-kjpaw_psl.1.0.0.UPF"])
            type_guess, a_guess, b_guess = guess_mesh_type(psp.mesh.r, psp.mesh.rab)
            @test type_guess == "log_1"
            @test a_guess ≈ psp.mesh.dx
            @test b_guess ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
        end

        @testset "[logarithmic type 1] He.pbe-hgh.UPF" begin
            psp = load_psp_file(upf2_filepaths["He.pbe-hgh.UPF"])
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

    @testset "fast_sphericalbesselj" begin
        for x in (rand(100) .* 100)
            @test sphericalbesselj(0, x) ≈ fast_sphericalbesselj0(x) atol = 1e-8
        end

        for l in 0:5
            for x in (rand(100) .* 100)
                @test sphericalbesselj(l, x) ≈ fast_sphericalbesselj(l, x) atol = 1e-8
            end
        end
        for l in 6:10
            @test_throws ErrorException("The case l = $l is not implemented") fast_sphericalbesselj(l, 1.0)
            @test fast_sphericalbesselj(l, 0.0) == 0.0
        end
    end

    @testset "trapezoidal reference" begin
        # Taken from Python Numerical Methods 21.03 "Trapezoidal Rule"
        a = 0
        b = π
        n = 11
        h = (b - a) / (n - 1)
        x = range(a, b, n)
        f = sin.(x)
        @test trapezoid(f, h) == 1.9835235375094546
        @test trapezoid(f, fill(h, n)) ≈ 1.9835235375094546
    end

    @testset "simpson reference" begin
        # Taken from Python Numerical Methods 21.04 "Simpson's Rule"
        a = 0
        b = π
        n = 11
        h = (b - a) / (n - 1)
        x = range(a, b, n)
        f = sin.(x)
        @test simpson(f, h) == 2.0001095173150043
        @test simpson(f, fill(h, n)) ≈ 2.0001095173150043
    end

    @testset "odd linear grid sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 100001
            x0 = 0.00
            dx = π / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            fx = sin.(x)
            @test integrator(fx, dx) ≈ 2.0
            @test integrator(fx, fill(dx, n)) ≈ 2.0
            @test integrator(fx, dx) ≈ integrator(fx, fill(dx, n))
        end
    end

    @testset "even linear grid sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 100000
            x0 = 0.00
            dx = π / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            fx = sin.(x)
            @test integrator(fx, dx) ≈ 2.0
            @test integrator(fx, fill(dx, n)) ≈ 2.0
            @test integrator(fx, dx) ≈ integrator(fx, fill(dx, n))
        end
    end

    @testset "odd log grid sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 1000001
            i = collect(1:n)
            b = π / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            fx = sin.(x)
            @test integrator(fx, dx) ≈ 2.0 atol=1e-5 rtol=1e-5
        end
    end

    @testset "even log grid sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 1000000
            i = collect(1:n)
            b = π / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            fx = sin.(x)
            @test integrator(fx, dx) ≈ 2.0 atol=1e-5 rtol=1e-5
        end
    end

    # Reference value for ∫[0,π] exp(-x^2) sin(x) dx 
    ref = -(√π * (-2erfi(1/2) + erfi(1/2 - im*π) + erfi(1/2 + im*π))) / (4 * exp(1/4))

    @testset "odd linear grid exp sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 100001
            x0 = 0.00
            dx = π / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            fx = @. exp(-(x^2)) * sin(x)
            @test integrator(fx, dx) ≈ real(ref)
            @test integrator(fx, fill(dx, n)) ≈ real(ref)
        end
    end

    @testset "even linear grid exp sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 100000
            x0 = 0.00
            dx = π / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            fx = @. exp(-(x^2)) * sin(x)
            @test integrator(fx, dx) ≈ real(ref)
            @test integrator(fx, fill(dx, n)) ≈ real(ref)
        end
    end

    @testset "odd log grid exp sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 1000001
            i = collect(1:n)
            b = π / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            fx = @. exp(-(x^2)) * sin(x)
            @test integrator(fx, dx) ≈ real(ref) atol=1e-5 rtol=1e-5
        end
    end

    @testset "even log grid exp sin" begin
        @testset for integrator in (trapezoid, simpson)
            n = 1000000
            i = collect(1:n)
            b = π / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            fx = @. exp(-(x^2)) * sin(x)
            @test integrator(fx, dx) ≈ real(ref) atol=1e-5 rtol=1e-5
        end
    end
end
