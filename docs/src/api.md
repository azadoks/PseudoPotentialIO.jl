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
AnalyticalPsP
HghPsP{T}
```

## Functions

```@docs
load_psp_file
load_psp
format
element
is_norm_conserving
is_ultrasoft
is_paw
formalism
has_spin_orbit
relativistic_treatment
has_nlcc
valence_charge
max_angular_momentum
n_projector_radials
n_pseudo_orbital_radials
local_potential_real
local_potential_fourier
projector_coupling
projector_real
projector_fourier
pseudo_energy_correction
core_charge_density_real
core_charge_density_fourier
valence_charge_density_real
valence_charge_density_fourier
pseudo_orbital_real
pseudo_orbital_fourier
```
