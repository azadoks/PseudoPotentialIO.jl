import PseudoPotentialIO: simpson, rectangle, trapezoid
import PseudoPotentialIO: abinit_corrected_trapezoid, qe_simpson, cp90_simpson
import PseudoPotentialIO: linear_mesh, logarithmic_mesh1, logarithmic_mesh2

@testset "Quadrature" begin

    NONUNIFORM_INTEGRATORS = [simpson, rectangle, trapezoid, qe_simpson, cp90_simpson]
    INTEGRATORS = [simpson, rectangle, trapezoid, abinit_corrected_trapezoid, qe_simpson,
                   cp90_simpson]

    @testset "Simpson vs. Python reference" begin
        # Taken from Python Numerical Methods 21.04 "Simpson's Rule"
        a = 0
        b = float(π)
        n = 11
        h = (b - a) / (n - 1)
        x = range(a, b, length=n)
        i_start = firstindex(x)
        i_stop = lastindex(x)
        f(i) = sin(x[i])
        @test simpson(f, i_start, i_stop, h) == 2.0001095173150043
        @test simpson(f, i_start, i_stop, fill(h, n)) ≈ 2.0001095173150043
    end

    @testset "Odd linear grid sin" begin
        @testset for integrator in INTEGRATORS
            n = 100001
            x0 = 0.00
            dx = float(π) / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ 2.0
            if integrator in NONUNIFORM_INTEGRATORS
                @test integrator(f, i_start, i_stop, fill(dx, n)) ≈ 2.0
                @test integrator(f, i_start, i_stop, dx) ≈ integrator(f, i_start, i_stop, fill(dx, n))
            end
        end
    end

    @testset "Even linear grid sin" begin
        @testset for integrator in INTEGRATORS
            n = 100000
            x0 = 0.00
            dx = float(π) / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ 2.0
            if integrator in NONUNIFORM_INTEGRATORS
                @test integrator(f, i_start, i_stop, fill(dx, n)) ≈ 2.0
                @test integrator(f, i_start, i_stop, dx) ≈ integrator(f, i_start, i_stop, fill(dx, n))
            end
        end
    end

    @testset "Odd log grid sin" begin
        @testset for integrator in NONUNIFORM_INTEGRATORS
            n = 1000001
            i = collect(1:n)
            b = float(π) / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ 2.0 atol=1e-5 rtol=1e-5
        end
    end

    @testset "Even log grid sin" begin
        @testset for integrator in NONUNIFORM_INTEGRATORS
            n = 1000000
            i = collect(1:n)
            b = float(π) / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ 2.0 atol=1e-5 rtol=1e-5
        end
    end

    # Reference value for ∫[0,π] exp(-x^2) sin(x) dx 
    ref = -(√π * (-2erfi(1/2) + erfi(1/2 - im*π) + erfi(1/2 + im*π))) / (4 * exp(1/4))

    @testset "Odd linear grid exp sin" begin
        @testset for integrator in INTEGRATORS
            n = 100001
            x0 = 0.00
            dx = float(π) / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = exp(-(x[i]^2)) * sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ real(ref)
            if integrator in NONUNIFORM_INTEGRATORS
                @test integrator(f, i_start, i_stop, fill(dx, n)) ≈ real(ref)
            end
        end
    end

    @testset "Even linear grid exp sin" begin
        @testset for integrator in INTEGRATORS
            n = 100000
            x0 = 0.00
            dx = float(π) / n
            x = linear_mesh.(1:n, dx, x0 - dx)
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = exp(-(x[i]^2)) * sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ real(ref)
            if integrator in NONUNIFORM_INTEGRATORS
                @test integrator(f, i_start, i_stop, fill(dx, n)) ≈ real(ref)
            end
        end
    end

    @testset "Odd log grid exp sin" begin
        @testset for integrator in NONUNIFORM_INTEGRATORS
            n = 1000001
            i = collect(1:n)
            b = float(π) / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = exp(-(x[i]^2)) * sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ real(ref) atol=1e-5 rtol=1e-5
        end
    end

    @testset "Even log grid exp sin" begin
        @testset for integrator in NONUNIFORM_INTEGRATORS
            n = 1000000
            i = collect(1:n)
            b = float(π) / n
            a = log(π / b) / (n - 1)
            x = logarithmic_mesh1.(1:n, a, b)
            dx = @. a * b * exp.(a * (i - 1))
            i_start = firstindex(x)
            i_stop = lastindex(x)
            f(i) = exp(-(x[i]^2)) * sin(x[i])
            @test integrator(f, i_start, i_stop, dx) ≈ real(ref) atol=1e-5 rtol=1e-5
        end
    end
end
