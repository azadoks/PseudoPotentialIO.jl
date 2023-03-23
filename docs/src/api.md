# API reference

## Index

```@index
Pages = ["api.md"]
```

## Types

```@docs
PsPFile
UpfFile
Psp8File
HghFile

AbstractPsP

NumericPsP{T}
NormConservingPsP{T}
UltrasoftPsP{T}
ProjectorAugmentedWavePsP{T}

AnalyticalPsP
HghPsP{T}

AbstractPsPQuantity

PsPChargeDensity
ValenceChargeDensity
CoreChargeDensity

PsPProjector
BetaProjector
ChiProjector

PsPCoupling
BetaCoupling
AugmentationFunction

PsPPotential
LocalPotential
AugmentationFunction

EvaluationSpace
RealSpace
FourierSpace
```

## Functions

```@docs
load_psp_file
load_psp
load_family_psp_files
load_family_psps
list_families
list_family_psps
show_family_summary
show_family_periodic_table
show_family_table

identifier
element
format
max_angular_momentum
angular_momenta
relativistic_treatment
formalism
valence_charge
atomic_charge
is_norm_conserving
is_ultrasoft
is_paw
has_spin_orbit
n_radials
n_angulars
has_quantity
get_quantity
cutoff_radius
psp_quantity_evaluator
psp_energy_correction

PseudoPotentialIO.hankel_transform
PseudoPotentialIO.build_interpolator_real
PseudoPotentialIO.simpson
PseudoPotentialIO.fast_sphericalbesselj
```
