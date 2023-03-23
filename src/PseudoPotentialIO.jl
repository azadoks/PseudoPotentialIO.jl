module PseudoPotentialIO
using Artifacts
using EzXML
using Interpolations
using LazyArtifacts
using LinearAlgebra
using OffsetArrays
using Polynomials
using Printf
using Statistics
using SHA
using PrettyTables
using OrderedCollections

using PeriodicTable: PeriodicTable
import Base.Broadcast.broadcastable
import Bessels: gamma, sphericalbesselj
import SpecialFunctions: erf

## DocStringExtensions Templates
# TODO they don't seem to be working at the moment ?
using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         $(METHODLIST)
                                         """

@template TYPES = """
                  $(TYPEDEF)
                  $(DOCSTRING)
                  $(TYPEDFIELDS)
                  """

## Quantity flags
export AbstractPsPQuantity
export PsPChargeDensity
export ValenceChargeDensity
export CoreChargeDensity
export PsPProjector
export BetaProjector
export ChiProjector
export PsPCoupling
export BetaCoupling
export PsPPotential
export LocalPotential
export AugmentationFunction
export EvaluationSpace
export RealSpace
export FourierSpace
include("psp/quantities.jl")

## File datastructures and interface
export PsPFile
export format
export element
export is_norm_conserving
export is_ultrasoft
export is_paw
export formalism
export has_spin_orbit
export relativistic_treatment
export has_core_density
export valence_charge
export max_angular_momentum
export n_radials
include("file/file.jl")

export UpfFile
include("file/upf.jl")
include("file/upf1.jl")
include("file/upf2.jl")

export Psp8File
include("file/psp8.jl")

export HghFile
include("file/hgh.jl")

## Pseudopotential datastructures and interface
export AbstractPsP
export identifier
export element
export max_angular_momentum
export n_angulars
export valence_charge
export atomic_charge
export is_norm_conserving
export is_ultrasoft
export is_paw
export has_spin_orbit
export has_quantity
export get_quantity
export cutoff_radius
export psp_quantity_evaluator
export psp_energy_correction
export angular_momenta
export relativistic_treatment
export formalism
include("psp/psp.jl")

export NumericPsP
include("psp/numeric.jl")

export NormConservingPsP
include("psp/norm_conserving.jl")

export UltrasoftPsP
include("psp/ultrasoft.jl")

export ProjectorAugmentedWavePsP
include("psp/paw.jl")

export AnalyticalPsP
include("psp/analytical.jl")

export HghPsP
include("psp/hgh.jl")

## Loading/listing/showing functions
export list_families
export list_family_psps
include("io/list.jl")

export load_family_psp_files
export load_family_psps
export load_psp_file
export load_psp
include("io/load.jl")

export show_family_summary
export show_family_table
export show_family_periodic_table
include("io/show.jl")

## Deprecated loaders
export load_upf
export load_psp8
include("deprecated/upf.jl")
include("deprecated/psp8.jl")

## Miscellaneous
include("common/hankel_transform.jl")
include("common/mesh.jl")
include("common/quadrature.jl")
include("common/spherical_bessel.jl")
include("common/interpolation.jl")
include("common/truncation.jl")
end
