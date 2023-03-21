abstract type AbstractPsPQuantity end

abstract type PsPChargeDensity <: AbstractPsPQuantity end
struct ValenceChargeDensity <: PsPChargeDensity end
struct CoreChargeDensity <: PsPChargeDensity end

abstract type PsPProjector <: AbstractPsPQuantity end
struct BetaProjector <: PsPProjector end
struct ChiProjector <: PsPProjector end

abstract type PsPCoupling <: AbstractPsPQuantity end
struct BetaCoupling <: PsPCoupling end

abstract type PsPPotential <: AbstractPsPQuantity end
struct LocalPotential <: PsPPotential end

struct AugmentationFunction <: AbstractPsPQuantity end

abstract type EvaluationSpace end
struct RealSpace <: EvaluationSpace end
struct FourierSpace <: EvaluationSpace end

# abstract type RelativisticTreatment end
# struct ScalarRelativistic <: RelativisticTreatment end
# struct FullyRelativistic <: RelativisticTreatment end
