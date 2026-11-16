using Test
using PseudoPotentialIO
using PseudoPotentialIO: is_norm_conserving, is_ultrasoft, is_paw, 
    has_spin_orbit, has_model_core_charge_density, identifier, 
    format, element, functional, valence_charge, formalism
using Aqua
using PseudoPotentialData
using PeriodicTable
using JSON
using LinearAlgebra
using Random

Random.seed!(0)

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

include("fixtures.jl")

@testset "PseudoPotentialIO.jl" begin
    if any(in.(("all", "aqua"), Ref(TAGS)))
        include("aqua.jl")
    end

    if any(in.(("all", "file"), Ref(TAGS)))
        include("file/file.jl")
    end

    if any(in.(("all", "file", "upf"), Ref(TAGS)))
        include("file/upf/upf.jl")
    end

    if any(in.(("all", "file", "upf", "upf1"), Ref(TAGS)))
        include("file/upf/upf1.jl")
        include("file/upf/file_dict_v1.jl")
    end

    if any(in.(("all", "file", "upf", "upf2"), Ref(TAGS)))
        include("file/upf/upf2.jl")
        include("file/upf/file_dict_v2.jl")
    end

    if any(in.(("all", "file", "psml"), Ref(TAGS)))
        include("file/psml/psml.jl")
        include("file/psml/file_dict.jl")
    end

    if any(in.(("all", "file", "psp8"), Ref(TAGS)))
        include("file/psp8/psp8.jl")
        include("file/psp8/file_dict.jl")
    end

    if any(in.(("all", "file", "hgh"), Ref(TAGS)))
        include("file/hgh/hgh.jl")
        include("file/hgh/file_dict.jl")
    end

    if any(in.(("all", "file", "validation"), Ref(TAGS)))
        include("file/validation.jl")
    end

    if any(in.(("all", "file", "roundtrip"), Ref(TAGS)))
        include("file/roundtrip.jl")
    end

    if any(in.(("all", "file", "upf", "threads"), Ref(TAGS)))
        include("file/threads.jl")
    end
end
