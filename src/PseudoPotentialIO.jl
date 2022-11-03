module PseudoPotentialIO

using EzXML
using DocStringExtensions
using Printf
using Interpolations
using Statistics
import PeriodicTable

include("common.jl")

export AbstractPsP
export format
export element
export l_max
export n_proj_radial
export n_proj
export n_pseudo_wfc
export z_valence
export is_paw
export is_ultrasoft
export is_norm_conserving
export is_coulomb
export formalism
export has_spin_orbit
export relativistic_treatment
export has_nlcc
export e_kb
export v_local_real
export v_local_correction_real
export v_local_correction_fourier
export v_local_fourier
export projector_real
export projector_fourier
export pseudo_energy_correction
include("interface.jl")

export UpfPsP
include("upf.jl")
include("upf1.jl")
include("upf2.jl")

export Psp8PsP
include("psp8.jl")

export HghPsP
include("hgh.jl")

end
