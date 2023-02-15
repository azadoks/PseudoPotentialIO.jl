using Test
using PseudoPotentialIO
using Aqua
using LazyArtifacts
using Random
using PeriodicTable
using SpecialFunctions
using QuadGK
using Quaternions
using Interpolations

Random.seed!(0)

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

include("fixtures.jl")

@testset "PseudoPotentialIO.jl" begin
    if any(in.(("all", "aqua"), Ref(TAGS)))
        include("aqua.jl")
    end

    if any(in.(("all", "common"), Ref(TAGS)))
        include("common.jl")
    end

    if any(in.(("all", "file"), Ref(TAGS)))
        include("file/file.jl")
    end

    if any(in.(("all", "file", "upf"), Ref(TAGS)))
        include("file/upf.jl")
    end

    if any(in.(("all", "file", "upf", "upf1"), Ref(TAGS)))
        include("file/upf1.jl")
    end

    if any(in.(("all", "file", "upf", "upf2"), Ref(TAGS)))
        include("file/upf2.jl")
    end

    if any(in.(("all", "file", "psp8"), Ref(TAGS)))
        include("file/psp8.jl")
    end

    if any(in.(("all", "file", "hgh"), Ref(TAGS)))
        include("file/hgh.jl")
    end

    if any(in.(("all", "psp"), Ref(TAGS)))
        include("psp/psp.jl")
    end

    if any(in.(("all", "psp", "numeric"), Ref(TAGS)))
        include("psp/numeric.jl")
    end

    if any(in.(("all", "psp", "numeric", "norm_conserving"), Ref(TAGS)))
        include("psp/norm_conserving.jl")
    end

    if any(in.(("all", "psp", "numeric", "ultrasoft"), Ref(TAGS)))
        include("psp/ultrasoft.jl")
    end

    if any(in.(("all", "psp", "analytical"), Ref(TAGS)))
        include("psp/analytical.jl")
    end

    if any(in.(("all", "psp", "analytical", "upf_hgh"), Ref(TAGS)))
        include("psp/upf_hgh.jl")
    end

    if any(in.(("all", "psp", "numeric", "upf_psp8"), Ref(TAGS)))
        include("psp/upf_psp8.jl")
    end
end
