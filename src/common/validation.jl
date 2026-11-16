#= Validation of the file structs against the specification of their format. The
dictionary readers (`file_dict.jl`) take what the file has; the constructors of the structs
check that it follows the specification (`validate`), and treat what does not as `VALIDATION`
says. =#

"""
What the constructors of the file structs (e.g. [`UpfFile`](@ref)) do with a value that does
not follow the specification of the format:

- `:error`: throw an `ArgumentError`;
- `:warn` (the default): log a warning, and keep the value;
- `:none`: keep the value.

Some checks are looser than the specification, where common files (or the reference
reader of the format) do not follow it: those warn, but never throw (`@soft_check`).

Set it for a block of code with [`with_validation`](@ref), or for one file with the
`validation` keyword of [`load_psp_file`](@ref) and of the constructors from a file (e.g.
`UpfFile(path; validation=:error)`). It is a `ScopedValue`: the tasks started in the block
see it too. What cannot be stored at all (text where a number is needed, a missing required
value, arrays whose sizes do not fit) is always an error.
"""
const VALIDATION = ScopedValue(:warn)

# The identifier of the file being read, for the messages
const VALIDATION_SOURCE = ScopedValue("")

const VALIDATION_MODES = (:error, :warn, :none)

"""
Run `f()` with the validation mode `mode` (see [`VALIDATION`](@ref)); the messages name
`source` (e.g. the file).
"""
function with_validation(f, mode::Symbol; source::AbstractString=VALIDATION_SOURCE[])
    mode in VALIDATION_MODES ||
        throw(ArgumentError("validation must be one of $VALIDATION_MODES, not :$mode"))
    return with(f, VALIDATION => mode, VALIDATION_SOURCE => String(source))
end

# `soft`: a warning, also in the mode `:error`
function spec_violation(format::AbstractString, message::AbstractString; soft=false)
    mode = VALIDATION[]
    mode === :none && return nothing
    source = VALIDATION_SOURCE[]
    text = isempty(source) ? "$format: $message" : "$format ($source): $message"
    mode === :warn || soft ? (@warn text) : throw(ArgumentError(text))
    return nothing
end

"""
    @check format condition message

A specification of `format` (e.g. "UPF"): unless `condition` holds, a violation with
`message` (only built then), which `VALIDATION` turns into an error, a warning or nothing.
"""
macro check(format, condition, message)
    return :($(esc(condition)) || spec_violation($(esc(format)), $(esc(message))))
end

"""
    @soft_check format condition message

As `@check`, for what the specification asks but files commonly break: a warning (unless
`VALIDATION` is `:none`), never an error.
"""
macro soft_check(format, condition, message)
    return :($(esc(condition)) || spec_violation($(esc(format)), $(esc(message)); soft=true))
end

# The checks of a struct `x` (with `@check`); nothing by default
validate(x) = nothing

# `x`, once validated: the inner constructors are `T(fields...) = validated(new(fields...))`
validated(x) = (validate(x); x)

# Spin-orbit: j = l ± 1/2 (and j > 0)
valid_j(l, j) = j > 0 && abs(j - l) == 0.5
