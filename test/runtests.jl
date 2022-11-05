using Test
using PseudoPotentialIO
using Random
using PeriodicTable
using SpecialFunctions
using Interpolations
using QuadGK

Random.seed!(0)

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

@testset "PseudoPotentialIO.jl" begin
    if "all" in TAGS || "common" in TAGS
        include("common.jl")
    end

    if "all" in TAGS || "upf" in TAGS || "numeric" in TAGS
        include("upf.jl")
    end

    if "all" in TAGS || "psp8" in TAGS || "numeric" in TAGS
        include("psp8.jl")
    end

    if "all" in TAGS || "hgh" in TAGS || "analytical" in TAGS
        include("hgh.jl")
    end
end
