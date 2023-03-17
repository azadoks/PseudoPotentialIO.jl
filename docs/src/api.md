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
```

## Functions

```@docs
load_psp_file
load_psp
list_families
show_family_table
show_family_list
list_psp

identifier
element
format
max_angular_momentum
angular_momenta
relativistic_treatment
formalism
valence_charge
atomic_charge
n_projector_radials
n_projector_angulars
n_chi_function_radials
n_chi_function_angulars
is_norm_conserving
is_ultrasoft
is_paw
has_spin_orbit
has_core_density
has_valence_density
has_chi_functions
local_potential_cutoff_radius
projector_cutoff_radius
chi_function_cutoff_radius
valence_charge_density_cutoff_radius
core_charge_density_cutoff_radius
projector_radial_indices
chi_function_radial_indices
projector_coupling
pseudo_energy_correction
augmentation_coupling
local_potential_real
local_potential_fourier
projector_real
projector_fourier
chi_function_real
chi_function_fourier
valence_charge_density_real
valence_charge_density_fourier
core_charge_density_real
core_charge_density_fourier
augmentation_real
augmentation_fourier

PseudoPotentialIO.hankel_transform
PseudoPotentialIO.build_interpolator_real
PseudoPotentialIO.simpson
PseudoPotentialIO.fast_sphericalbesselj
```
