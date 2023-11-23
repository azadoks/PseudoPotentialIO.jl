abstract type AbstractPsPQuantity end

abstract type AbstractProjector <: AbstractPsPQuantity end
struct NumericProjector{T,A<:AbstractVector{T}} <: AbstractProjector
    n::Int
    l::Int
    j::T
    r::RadialMesh{T}
    r²f::A
    sphericalbesselj::Function
end

abstract type AbstractState <: AbstractPsPQuantity end
struct NumericState{T,A<:AbstractVector{T}} <: AbstractState
    n::Int
    l::Int
    j::T
    r::RadialMesh{T}
    r²f::A
    sphericalbesselj::Function
end

abstract type AbstractNumericLocalPotential <: AbstractPsPQuantity end
struct NumericLocalPotential{T,A<:AbstractVector{T}} <: AbstractPsPQuantity
    r::RadialMesh{T}
    f::A
end

abstract type AbstractDensity <: AbstractPsPQuantity end
struct NumericDensity{T,A<:AbstractVector{T}} <: AbstractPsPQuantity
    r::RadialMesh{T}
    r²f::A
end

abstract type AbstractAugmentation <: AbstractPsPQuantity end
struct NumericAugmentation{T,A<:AbstractVector{T}} <: AbstractPsPQuantity
    n::Int
    m::Int
    l::Int
    r::RadialMesh{T}
    r²f::A
    sphericalbesselj::Function
end

# # # TODO: add l, n, ... indices as components of the flagging structure
# # # TODO: this should greatly reduce code duplication!

# """
# Abstact type used for dispatch on different pseudopotential quantities"
# """
# abstract type AbstractPsPQuantity end

# """
# Generic atomic charge denstiy quantity
# """
# abstract type AbstractDensity <: AbstractPsPQuantity end
# """
# Pseudo-atomic valence charge density
# """
# struct ValenceChargeDensity <: AbstractDensity end
# """
# Pseudo-atomic core charge density
# """
# struct CoreChargeDensity <: AbstractDensity end

# """
# Generic projector/orbital-like quantity
# """
# abstract type AbstractProjector <: AbstractPsPQuantity end
# @doc raw"""
# Kleinman-Bylander non-local ``\beta`` projector
# """
# struct NumericProjector <: AbstractProjector end
# @doc raw"""
# ``\chi`` projector
# """
# struct NumericState <: AbstractProjector end

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

# """
# Generic potential quantity
# """
# abstract type AbstractPotential <: AbstractPsPQuantity end
# """
# Local potential
# """
# struct NumericLocalPotential <: AbstractPotential end

abstract type NumericLocalPotentialCorrection <: AbstractPsPQuantity end
struct ErfCorrection <: NumericLocalPotentialCorrection end
struct CoulombCorrection <: NumericLocalPotentialCorrection end

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
