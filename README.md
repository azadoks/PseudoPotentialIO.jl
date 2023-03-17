# PseudoPotentialIO.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://azadoks.github.io/PseudoPotentialIO.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://azadoks.github.io/PseudoPotentialIO.jl/dev)
[![Build Status](https://github.com/azadoks/PseudoPotentialIO.jl/workflows/CI/badge.svg)](https://github.com/azadoks/PseudoPotentialIO.jl/actions)

PseudoPotentialIO aims to provide parsers for common pseudopotential file formats used by density functional theory codes and an interface for accesssing the quantities that they contain.

The following file formats are supported or have planned (or no planned) support.
If your favorite format does not appear in the table below, please file an [issue](https://github.com/azadoks/PseudoPotentialIO.jl/issues)!

| Format                   | Read        | Standardize                   |
|--------------------------|-------------|-------------------------------|
| UPF (old)                | ✅          | `NormConserving`, `Ultrasoft` |
| UPF (v2.0.1 with schema) | ✅          | `NormConserving`, `Ultrasoft` |
| PSP1                     | Not planned |                               |
| PSP3/HGH                 | Planned     | `Hgh`                         |
| PSP4/Teter               | Not planned |                               |
| PSP5/"Phoney"            | Not planned |                               |
| PSP6/FHI98PP             | Not planned |                               |
| PSP7/ABINIT PAW          | Not planned |                               |
| PSP8/ONCVPSP             | ✅          | `NormConserving`              |
| PSP9/PSML                | Planned     |                               |
| PSP17/ABINIT PAW XML     | Planned     |                               |
| FPMD XML                 | Not planned |                               |
| Vanderbilt USPP          | Planned     |                               |

Pseudopotentials are read into `[Format]File` structs which mirror very closely the contents of the file.
However, different file formats provide important physical quantities in slightly different forms.
For a more consistent interface to the physical quantities with consitent units, meanings, and dimensions, `*File` structs can be converted to `*Psp` pseudopotential structs.

The following type tree for representing different types of pseudopotentials is implemented

```
AbstractPsP
|
--- NumericPsP (pseudopotentials on grids)
|   |
|   --- NormConservingPsP (UPF, PSP8)
|   |
|   --- UltrasoftPsP (UPF)
|   |
|   --- [Placeholder] ProjectorAugmentedWavePsP (UPF)
|
--- AnalyticPsP (pseudopotentials with analytic forms)
    |
    --- HghPsP (HGH)
```

Unsupported features:
- Fully relativistic (spin-orbit coupling) pseudopotentials
- Semilocal pseudopotentials

Supported features:
- Non-linear core corrections
- Multiprojector pseudopotentials
- Linear and logarithmic meshes
