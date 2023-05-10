# TODO: add l, n, ... indices as components of the flagging structure
# TODO: this should greatly reduce code duplication!

"""
Abstact type used for dispatch on different pseudopotential quantities"
"""
abstract type AbstractPsPQuantity end

"""
Generic atomic charge denstiy quantity
"""
abstract type PsPChargeDensity <: AbstractPsPQuantity end
"""
Pseudo-atomic valence charge density
"""
struct ValenceChargeDensity <: PsPChargeDensity end
"""
Pseudo-atomic core charge density
"""
struct CoreChargeDensity <: PsPChargeDensity end

"""
Generic projector/orbital-like quantity
"""
abstract type PsPProjector <: AbstractPsPQuantity end
@doc raw"""
Kleinman-Bylander non-local ``\beta`` projector
"""
struct BetaProjector <: PsPProjector end
@doc raw"""
``\chi`` projector
"""
struct ChiProjector <: PsPProjector end

"""
Generic coupling (matrix) quantity
"""
abstract type PsPCoupling <: AbstractPsPQuantity end
"""
Kleinman-Bylander energies / ``D`` coupling matrix
"""
struct BetaCoupling <: PsPCoupling end
@doc raw"""
Ultrasoft/PAW augmentation function coupling matrix
"""
struct AugmentationCoupling <: PsPCoupling end

"""
Generic potential quantity
"""
abstract type PsPPotential <: AbstractPsPQuantity end
"""
Local potential
"""
struct LocalPotential <: PsPPotential end

abstract type LocalPotentialCorrection <: AbstractPsPQuantity end
struct ErfCorrection <: LocalPotentialCorrection end
struct CoulombCorrection <: LocalPotentialCorrection end

@doc raw"""
Ultrasoft / PAW augmentation function ``q`` or ``Q``
"""
struct AugmentationFunction <: AbstractPsPQuantity end

"""
Abstact type used in dispatch selection of evaluation space
"""
abstract type EvaluationSpace end
"""
Evaluate a quantity in real space
"""
struct RealSpace <: EvaluationSpace end
"""
Evaluate the Fourier transform of a quantity, i.e. evaluate it in Fourier space
"""
struct FourierSpace <: EvaluationSpace end

# abstract type RelativisticTreatment end
# struct ScalarRelativistic <: RelativisticTreatment end
# struct FullyRelativistic <: RelativisticTreatment end
