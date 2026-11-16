# PseudoPotentialIO.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/PseudoPotentialIO.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/PseudoPotentialIO.jl/dev)
[![Build Status](https://github.com/JuliaMolSim/PseudoPotentialIO.jl/workflows/CI/badge.svg)](https://github.com/JuliaMolSim/PseudoPotentialIO.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaMolSim/PseudoPotentialIO.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaMolSim/PseudoPotentialIO.jl)

PseudoPotentialIO aims to provide parsers for common pseudopotential file formats used by density functional theory codes and an interface for accesssing the quantities that they contain.

## Reading pseudopotentials

The following file formats are supported (✅) or not (-).
If your favorite format does not appear in the table below, please file an [issue](https://github.com/JuliaMolSim/PseudoPotentialIO.jl/issues)!

| Format                        | Read        | Write |
|-------------------------------|-------------|-------|
| UPF v1                        | ✅          | -     |
| UPF v2.0.1                    | ✅          | ✅    |
| UPF (`qe_pp` schema)          | -           | -     |
| PSP1/Teter                    | -           | -     |
| PSP2/GTH (1996)               | -           | -     |
| PSP3/HGH                      | ✅          | ✅    |
| PSP4/Teter                    | -           | -     |
| PSP5/"Phoney"                 | -           | -     |
| PSP6/FHI98PP                  | -           | -     |
| PSP7/ABINIT PAW               | -           | -     |
| PSP8/ONCVPSP                  | ✅          | ✅    |
| PSP9/PSML                     | ✅          | ✅    |
| PSP10/HGH (full h/k)          | ✅          | ✅    |
| PSP17/PAW XML (ABINIT, GPAW)  | -           | -     |
| GTH/HGH (CP2K)                | ✅          | ✅    |
| FPMD XML                      | -           | -     |
| SIESTA `.psf`/`.vps`          | -           | -     |
| CASTEP `.usp`/`.recpot`       | -           | -     |
| CPMD `.psp`                   | -           | -     |
| Vanderbilt USPP               | -           | -     |
| VASP POTCAR                   | -           | -     |

Pseudopotentials are read by calling `load_psp_file(path)`.
This returns one of the `[Format]File` structs which mirror very closely the contents of the file.
However, different file formats provide important physical quantities in slightly different forms.
Be sure to check the docstrings and comments here in PseudoPotentialIO.jl or the relevant documentation of the file format.

The structs check that the file follows the specification of its format (for example, that the
flags of a UPF header agree with the data in the file).
By default, what does not is a warning; make it an error with
`load_psp_file(path; validation=:error)`, or turn the checks off with `validation=:none`.
Where common files break the specification (e.g. QE's own files, ONCVPSP's PSML), the checks
only warn, also with `:error`.
`PseudoPotentialIO.with_validation(f, mode)` sets the mode for a whole block of code.

## Writing pseudopotentials
This package also supports writing pseudopotentials, currently only in
the UPF v2.0.1, PSP8, HGH and PSML (1.1) formats. This is done using `save_psp_file`:
```jl
save_psp_file("path/to/pseudo.upf", pseudo)
```
where `pseudo` is a `UpfFile`, `Psp8File`, `HghFile` or `PsmlFile`.
Support for additional formats might be added in the future.

## Converting pseudopotentials
Conversion from PSP8 to UPF is implemented as well,
by calling the `UpfFile` constructor with a `Psp8File`.
Conversions between other formats might be added in the future.

Here is a full example of reading a PSP8 file,
converting it to a UPF file and then writing it:
```jl
psp8 = load_psp_file("path/to/pseudo.psp8")
save_psp_file("path/to/pseudo.upf", UpfFile(psp8))
```
