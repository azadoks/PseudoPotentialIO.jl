import PseudoPotentialIO: Trapezoid, Simpson, AbinitCorrectedTrapezoid, QESimpson
import PseudoPotentialIO: hankel_transform
# TODO: test hankel
@testset "Hankel transform" begin
    quadrature_methods = [Trapezoid(), Simpson(), AbinitCorrectedTrapezoid(), QESimpson()]
    r²f = nothing
    r = 0:0.1:10
    dr = 0.1
    l = 0
    for quadrature_method in quadrature_methods, l in 0:5
        @test isnothing(hankel_transform(r, r²f, dr, l; quadrature_method)(0.0))
    end
end
