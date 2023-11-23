# TODO: test truncation
import PseudoPotentialIO: find_truncation_index
@testset "Truncation index" begin
    @test isnothing(find_truncation_index(nothing, nothing))
    @test isnothing(find_truncation_index(nothing, 1e-2))
    
    @test find_truncation_index([], nothing) == 0
    @test find_truncation_index([100, 9, 1], nothing) == 3    

    @test find_truncation_index([100], 10) == 1
    @test find_truncation_index([1], 10) == 1
    
    @test find_truncation_index([100, 11], 10) == 2
    @test find_truncation_index([100, 9], 10) == 1
    
    @test find_truncation_index([100, 11, 1], 10) == 2
    @test find_truncation_index([100, 9, 1], 10) == 1
    @test find_truncation_index([100, 11, 9, 1], 10) == 2
    @test find_truncation_index([100, 10, 8, 1], 10) == 1

    @test find_truncation_index([100, 50, 0, 50, 10, 9, 8], 10) == 4
end
