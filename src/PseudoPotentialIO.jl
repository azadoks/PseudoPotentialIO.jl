module PseudoPotentialIO
using Artifacts
using DocStringExtensions
using EzXML
using Interpolations
using LazyArtifacts
using LinearAlgebra
using OffsetArrays
using Polynomials
using Printf
using SpecialFunctions
using Statistics

import Base.Broadcast.broadcastable
import PeriodicTable

## File I/O
export PsPFile
export format
export element
export is_norm_conserving
export is_ultrasoft
export is_paw
export formalism
export has_spin_orbit
export relativistic_treatment
export has_nlcc
export valence_charge
export max_angular_momentum
export n_projectors
export n_pseudo_orbitals
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
export element
export formalism
export relativistic_treatment
export has_spin_orbit
export has_nlcc
export valence_charge
export max_angular_momentum
export n_projectors
export n_pseudo_orbitals
export local_potential_real
export local_potential_fourier
export projector_coupling
export projector_real
export projector_fourier
export pseudo_energy_correction
export core_charge_density_real
export core_charge_density_fourier
export valence_charge_density_real
export valence_charge_density_fourier
export pseudo_orbital_real
export pseudo_orbital_fourier
include("psp/psp.jl")

export NumericPsP
include("psp/numeric.jl")

export NormConservingPsP
include("psp/norm_conserving.jl")

export UltrasoftPsP
include("psp/ultrasoft.jl")

export AnalyticalPsP
include("psp/analytical.jl")

export HghPsP
include("psp/hgh.jl")

## Core functions
export load_psp_file
export load_psp
include("load.jl")

## Miscellaneous
include("common/mesh.jl")
include("common/quadrature.jl")
include("common/spherical_bessel.jl")
end
