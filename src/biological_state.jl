# Prototype for OceanBioME: transform the state seen by biology without
# changing prognostic tracer fields.

abstract type AbstractBiologicalState end

"""Use prognostic tracer values unchanged when evaluating biogeochemical tendencies."""
struct IdentityBiologicalState <: AbstractBiologicalState end

"""Present nonnegative tracer values to biology while leaving prognostic fields unchanged."""
struct NonnegativeBiologicalState <: AbstractBiologicalState end

# Trait for signed state variables used by biology. Biology-owned concentrations
# are nonnegative by default; models only override variables that legitimately
# cross zero (for example temperature).
@inline biologically_nonnegative_tracer(::Type, ::Val) = Val(true)

@inline biological_value(::IdentityBiologicalState, model_type, val_name, value) = value
@inline biological_value(::NonnegativeBiologicalState, model_type, val_name, value) =
    biological_value(biologically_nonnegative_tracer(model_type, val_name), value)

@inline biological_value(::Val{true}, value) = max(value, zero(value))
@inline biological_value(::Val{false}, value) = value

# Continuous-form models receive scalar tracer values followed, in some call
# paths, by auxiliary values. Transform only the required tracer prefix and
# pass any remaining values through unchanged.
@inline biological_values(state, model_type, ::Tuple{}, values::Tuple) = values

@inline function biological_values(state, model_type, names::Tuple, values::Tuple)
    name = first(names)
    value = first(values)
    transformed = biological_value(state, model_type, Val(name), value)
    return (transformed, biological_values(state, model_type, Base.tail(names), Base.tail(values))...)
end

@inline function biological_values(state, bgc, values::Tuple)
    names = required_biogeochemical_tracers(bgc)
    return biological_values(state, typeof(bgc), names, values)
end

# Discrete-form models receive Field objects and index them inside their
# tendency functions. `NonnegativeBiologicalFields` is a lazy proxy, so only
# fields actually read by a tendency get wrapped.
struct NonnegativeField{F}
    field :: F
end

@inline Base.getindex(field::NonnegativeField, I...) = begin
    value = @inbounds getfield(field, :field)[I...]
    max(value, zero(value))
end

Base.eltype(field::NonnegativeField) = eltype(getfield(field, :field))
Base.eltype(::Type{NonnegativeField{F}}) where F = eltype(F)
Base.size(field::NonnegativeField, args...) = size(getfield(field, :field), args...)
Base.axes(field::NonnegativeField, args...) = axes(getfield(field, :field), args...)
Base.parent(field::NonnegativeField) = getfield(field, :field)

# Forward ordinary field properties such as `grid`, `data`, etc. This does not
# change dispatch for functions explicitly requiring an Oceananigans Field;
# that is one of the items the prototype tests need to validate.
@inline function Base.getproperty(field::NonnegativeField, name::Symbol)
    name === :field && return getfield(field, :field)
    return getproperty(getfield(field, :field), name)
end

@inline biological_field(::Val{true}, field) = NonnegativeField(field)
@inline biological_field(::Val{false}, field) = field

struct NonnegativeBiologicalFields{B, F}
    fields :: F
end

NonnegativeBiologicalFields(bgc::B, fields::F) where {B, F} =
    NonnegativeBiologicalFields{B, F}(fields)

Base.keys(fields::NonnegativeBiologicalFields) = keys(getfield(fields, :fields))
Base.propertynames(fields::NonnegativeBiologicalFields) = propertynames(getfield(fields, :fields))

@inline function Base.getproperty(fields::NonnegativeBiologicalFields{B}, name::Symbol) where B
    name === :fields && return getfield(fields, :fields)
    field = getproperty(getfield(fields, :fields), name)
    return biological_field(biologically_nonnegative_tracer(B, Val(name)), field)
end

@inline Base.getindex(fields::NonnegativeBiologicalFields, name::Symbol) = getproperty(fields, name)

@inline biological_fields(::IdentityBiologicalState, bgc, fields) = fields
@inline biological_fields(::NonnegativeBiologicalState, bgc, fields) =
    NonnegativeBiologicalFields(bgc, fields)
