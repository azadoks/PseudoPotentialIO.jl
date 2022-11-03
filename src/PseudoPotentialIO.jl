module PseudoPotentialIO

using EzXML
using DocStringExtensions
using Printf
using Interpolations
import PeriodicTable

include("common.jl")

export AbstractPsP
include("interface.jl")

export UpfPsP
include("upf.jl")
include("upf1.jl")
include("upf2.jl")

export Psp8PsP
include("psp8.jl")

end
