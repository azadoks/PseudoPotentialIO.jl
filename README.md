# PseudoPotentialIO.jl

PseudoPotentialIO aims to provide parsers for common pseudopotential file formats used by density functional theory codes and an interface for accesssing the quantities that they contain.

The following file formats are supported or have planned (or no planned) support.
If your favorite format does not appear in the table below, please file an [issue](https://github.com/azadoks/PseudoPotentialIO.jl/issues)!

|                          | Read |
|--------------------------|------|
| UPF (old)                | ✅   |
| UPF (v2.0.1 with schema) | ✅   |
| PSP1                     | Not planned     |
| PSP3/HGH                 | Planned     |
| PSP4/Teter               | Not planned     |
| PSP5/"Phoney"            | Not planned     |
| PSP6/FHI98PP             | Not planned     |
| PSP8/ONCVPSP             | ✅   |
| PSML                     | Planned     |
| ATOMPAW                  | Planned     |

<!-- The pseudopotential interface provides accessors for the following quantities:

- Valence charge
- Element
- Dimensions (mesh size, number of non-local projectors, maximum angular momentum, etc.)
- Local part of the potential in real and Fourier space
- Non-local projectors in real and Fourier space
- Kleinman-Bylander energies
- Model core charge density in real and Fourier space (for non-linear core correction)
- Valence charge density in real and Fourier space
- Pseudo-atomic orbitals in real and Fourier space
 -->
