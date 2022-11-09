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
    if any(in.(("all", "common"), Ref(TAGS)))
        include("common.jl")
    end

    if any(in.(("all", "file", "upf"), Ref(TAGS)))
        include("file/upf.jl")
    end

    if any(in.(("all", "file", "psp8"), Ref(TAGS)))
        include("file/psp8.jl")
    end

    if any(in.(("all", "file", "hgh"), Ref(TAGS)))
        include("file/hgh.jl")
    end

    if any(in.(("all", "psp", "numeric", "norm_conserving"), Ref(TAGS)))
        include("psp/norm_conserving.jl")
    end

    if any(in.(("all", "psp", "numeric", "ultrasoft"), Ref(TAGS)))
        include("psp/ultrasoft.jl")
    end

    if any(in.(("all", "psp", "analytical", "hgh"), Ref(TAGS)))
        include("psp/hgh.jl")
    end
end
