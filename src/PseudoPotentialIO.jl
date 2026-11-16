module PseudoPotentialIO
using OrderedCollections: OrderedDict
using XML: XML
using Printf
using PeriodicTable
using Statistics
using ScopedValues: ScopedValue, with

using PeriodicTable: PeriodicTable
import Base.Broadcast.broadcastable

## DocStringExtensions Templates
using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         $(METHODLIST)
                                         """

@template TYPES = """
                  $(TYPEDEF)
                  $(DOCSTRING)
                  $(TYPEDFIELDS)
                  """
## Data
include("data/libxc_functionals.jl")
include("data/psp8_functionals.jl")
include("data/upf_functionals.jl")

## Common utilities
include("common/mesh.jl")
include("common/fortran.jl")
include("common/xml.jl")
include("common/validation.jl")

## File datastructures and interface
# These are only 'public' to avoid clashes in the ecosystem with functions
# of similar names, e.g. in Libxc.jl and AtomsBase.jl
export PsPFile
#= TODO: Can be changed to this once 1.11 is hard enforced
public is_norm_conserving
public is_ultrasoft
public is_paw
public has_spin_orbit
public has_model_core_charge_density
public identifier
public format
public element
public functional
public valence_charge
public formalism
public with_validation
=#
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public identifier, format, element, functional, valence_charge, formalism, is_norm_conserving, is_ultrasoft, is_paw, has_spin_orbit, has_model_core_charge_density, with_validation, VALIDATION"))
end
include("file/file.jl")

# Each format: `file_dict.jl` reads and writes the file as nested dictionaries,
# `struct.jl` defines its `PsPFile`, `dict_struct.jl` converts between the two
export UpfFile
include("file/upf/struct.jl")
include("file/upf/file_dict_v1.jl")  # read into the dictionaries of UPF v2
include("file/upf/file_dict_v2.jl")
include("file/upf/dict_struct.jl")

export PsmlFile
include("file/psml/struct.jl")
include("file/psml/file_dict.jl")
include("file/psml/dict_struct.jl")

export Psp8File
include("file/psp8/struct.jl")
include("file/psp8/file_dict.jl")
include("file/psp8/dict_struct.jl")

export HghFile
include("file/hgh/struct.jl")
include("file/hgh/file_dict_cp2k.jl")
include("file/hgh/file_dict_abinit.jl")
include("file/hgh/dict_struct.jl")

## Loading/listing functions
export load_psp_file
include("load.jl")

## Save to file
export save_psp_file
include("save.jl")

# Conversion
include("conversion/to_upf.jl")

end
