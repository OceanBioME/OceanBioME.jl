using Oceananigans: fields, Simulation
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device, architecture, on_architecture
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry, required_biogeochemical_tracers
using Oceananigans.ImmersedBoundaries: immersed_cell

import Adapt: adapt_structure, adapt
import Base: summary, show
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!
import KernelAbstractions as KA

"""
    ClipNegativeTracers(; exclude = ())

Construct a negative-tracer treatment that clips negative tracer values to zero, excluding those listed in `exclude`.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
struct ClipNegativeTracers{E} end

ClipNegativeTracers(; exclude = ()) = ClipNegativeTracers{Tuple(exclude)}()

@inline excluded_tracers(::ClipNegativeTracers{E}) where E = E

function update_biogeochemical_state!(model, clip::ClipNegativeTracers)
    exclude = excluded_tracers(clip)
    for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

"""
    IgnoreNegativeTracerValues()

Construct a negative-tracer treatment that presents negative concentration-tracer values as zero
when evaluating biogeochemical processes while leaving prognostic tracer fields unchanged. Signed
environmental tracers such as temperature (`T`) and salinity (`S`) are evaluated unchanged.
"""
struct IgnoreNegativeTracerValues{N} end

IgnoreNegativeTracerValues() = IgnoreNegativeTracerValues{nothing}()

@inline resolve_negative_tracers(::IgnoreNegativeTracerValues{nothing}, bgc) =
    IgnoreNegativeTracerValues{required_biogeochemical_tracers(bgc)}()

update_biogeochemical_state!(model, ::IgnoreNegativeTracerValues) = nothing

# Evaluation-time treatments may be combined with state-mutating treatments in a
# tuple. Detect IgnoreNegativeTracerValues recursively so composition preserves
# its view-only semantics, including inside nested tuples produced while resolving
# model-dependent ScaleNegativeTracers treatments.
@inline ignores_negative_tracer_values(::IgnoreNegativeTracerValues) = Val(true)
@inline ignores_negative_tracer_values(::Any) = Val(false)
@inline ignores_negative_tracer_values(::Tuple{}) = Val(false)

@inline function ignores_negative_tracer_values(treatments::Tuple)
    first_ignored = ignores_negative_tracer_values(first(treatments))
    return ignores_negative_tracer_values(first_ignored, Base.tail(treatments))
end

@inline ignores_negative_tracer_values(::Val{true}, ::Tuple) = Val(true)
@inline ignores_negative_tracer_values(::Val{false}, treatments::Tuple) = ignores_negative_tracer_values(treatments)

# Concentration tracers are evaluated as nonnegative by default. Temperature and
# salinity are signed environmental state variables even when a biogeochemical
# model lists them among its required tracers. Models can add further exceptions.
@inline ignore_negative_value_for_tracer(::Type, ::Val) = Val(true)
@inline ignore_negative_value_for_tracer(::Type, ::Val{:T}) = Val(false)
@inline ignore_negative_value_for_tracer(::Type, ::Val{:S}) = Val(false)

# Continuous-form models receive required tracer values first, followed by any
# auxiliary values. The required tracer names are resolved into the treatment
# type before GPU execution so no Symbol-dependent dispatch is needed in kernels.
@generated function ignored_negative_values(::Type{B}, ::Val{names}, values::V) where {B, names, V <: Tuple}
    transformed = Any[]
    nrequired = length(names)
    nvalues = length(V.parameters)

    for i in 1:nvalues
        if i <= nrequired
            name = names[i]
            ignore = ignore_negative_value_for_tracer(B, Val(name))
            if ignore isa Val{true}
                push!(transformed, :(max(getfield(values, $i), zero(getfield(values, $i)))))
            else
                push!(transformed, :(getfield(values, $i)))
            end
        else
            push!(transformed, :(getfield(values, $i)))
        end
    end

    return :(($(transformed...),))
end

@inline ignore_negative_tracer_names(::IgnoreNegativeTracerValues{N}) where N = Val{N}()
@inline ignore_negative_tracer_names(::Any) = Val{nothing}()
@inline ignore_negative_tracer_names(::Tuple{}) = Val{nothing}()

@inline function ignore_negative_tracer_names(treatments::Tuple)
    first_names = ignore_negative_tracer_names(first(treatments))
    return ignore_negative_tracer_names(first_names, Base.tail(treatments))
end

@inline ignore_negative_tracer_names(names::Val{N}, ::Tuple) where N = names
@inline ignore_negative_tracer_names(::Val{nothing}, treatments::Tuple) = ignore_negative_tracer_names(treatments)

@inline evaluate_continuous_biogeochemistry(treatment, bgc, args...) =
    evaluate_continuous_biogeochemistry(ignore_negative_tracer_names(treatment), bgc, args...)

@inline function evaluate_continuous_biogeochemistry(names::Val{N}, bgc, val_name, x, y, z, t, values...) where N
    values = ignored_negative_values(typeof(bgc), names, values)
    return bgc(val_name, x, y, z, t, values...)
end

@inline evaluate_continuous_biogeochemistry(::Val{nothing}, bgc, args...) = bgc(args...)

# Discrete-form models index fields inside their tendency functions. Construct a
# NamedTuple with the same field names but nonnegative views for required
# concentration tracers. Field names are available as type parameters, so the
# transformation is fully static and GPU-safe.
struct NonnegativeValueField{F}
    field :: F
end

@inline function Base.getindex(field::NonnegativeValueField, I...)
    value = @inbounds getfield(field, :field)[I...]
    return max(value, zero(value))
end

Base.eltype(field::NonnegativeValueField) = eltype(getfield(field, :field))
Base.eltype(::Type{NonnegativeValueField{F}}) where F = eltype(F)
Base.size(field::NonnegativeValueField, args...) = size(getfield(field, :field), args...)
Base.axes(field::NonnegativeValueField, args...) = axes(getfield(field, :field), args...)
Base.parent(field::NonnegativeValueField) = getfield(field, :field)

@inline function Base.getproperty(field::NonnegativeValueField, name::Symbol)
    name === :field && return getfield(field, :field)
    return getproperty(getfield(field, :field), name)
end

# Chlorophyll-dependent light attenuation is evaluated through the same policy.
# Wrapping the derived chlorophyll field keeps this generic across light models
# while leaving the prognostic tracer fields untouched.
@inline chlorophyll(treatment, bgc, model) =
    chlorophyll(ignores_negative_tracer_values(treatment), bgc, model)
@inline chlorophyll(::Val{true}, bgc, model) = NonnegativeValueField(chlorophyll(bgc, model))
@inline chlorophyll(::Val{false}, bgc, model) = chlorophyll(bgc, model)

@generated function ignored_negative_fields(::Type{B}, ::Val{required}, fields::NamedTuple{names}) where {B, required, names}
    transformed = Any[]

    for (i, name) in enumerate(names)
        ignore = name in required && ignore_negative_value_for_tracer(B, Val(name)) isa Val{true}
        if ignore
            push!(transformed, :(NonnegativeValueField(getfield(fields, $i))))
        else
            push!(transformed, :(getfield(fields, $i)))
        end
    end

    return :(NamedTuple{$(QuoteNode(names))}(($(transformed...),)))
end

@inline evaluate_discrete_biogeochemistry(treatment, bgc, i, j, k, grid, val_name, clock, fields, auxiliary_fields) =
    evaluate_discrete_biogeochemistry(ignore_negative_tracer_names(treatment),
                                      bgc, i, j, k, grid, val_name, clock, fields, auxiliary_fields)

@inline function evaluate_discrete_biogeochemistry(names::Val{N}, bgc, i, j, k, grid, val_name, clock, fields, auxiliary_fields) where N
    fields = ignored_negative_fields(typeof(bgc), names, fields)
    return bgc(i, j, k, grid, val_name, clock, fields, auxiliary_fields)
end

@inline evaluate_discrete_biogeochemistry(::Val{nothing}, bgc, i, j, k, grid, val_name, clock, fields, auxiliary_fields) =
    bgc(i, j, k, grid, val_name, clock, fields, auxiliary_fields)

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{T, SA, FV, W}
      scalefactors :: SA
invalid_fill_value :: FV
              warn :: W

    function ScaleNegativeTracers(tracers, scalefactors::SA, invalid_fill_value::FV, warn::W) where {SA, W, FV}
        warn && error("Warning not currently implemented")
        return new{tracers, SA, FV, W}(scalefactors, invalid_fill_value, warn)
    end
end

@inline scaled_tracers(::ScaleNegativeTracers{T}) where T = T

@inline function Base.getproperty(scale::ScaleNegativeTracers, name::Symbol)
    name === :tracers && return scaled_tracers(scale)
    return getfield(scale, name)
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(scaled_tracers(snt),
                                                                      adapt(to, snt.scalefactors),
                                                                      adapt(to, snt.invalid_fill_value),
                                                                      adapt(to, snt.warn))

"""
    ScaleNegativeTracers(; invalid_fill_value = NaN, warn = false)

Construct a scaling treatment whose conserved tracer groups are inferred from the underlying
biogeochemical model when it is passed as `negative_tracers` to [`Biogeochemistry`](@ref) or a
model constructor.
"""
ScaleNegativeTracers(; invalid_fill_value = NaN, warn = false) =
    ScaleNegativeTracers(nothing, nothing, invalid_fill_value, warn)

@inline resolve_negative_tracers(scale::ScaleNegativeTracers{nothing}, bgc) =
    ScaleNegativeTracers(bgc; invalid_fill_value = scale.invalid_fill_value, warn = scale.warn)

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false, invalid_fill_value = NaN)

Constructs a negative-tracer treatment to scale `tracers` so that none are negative. Use like:
```julia
negative_tracers = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; negative_tracers)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`ClipNegativeTracers`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

Future plans include implement a positivity-preserving timestepping scheme as the ideal alternative.

~~If `warn` is true then scaling will raise a warning.~~

`invalid_fill_value` specifies the value to set the total cell content to if the total is less than 0
(meaning that total tracer conservation can not be enforced). If the value is set to anything other than
`NaN` this scheme no longer conserves mass. While this may be useful to prevent spurious numerics leading
to crashing care should be taken that the mass doesn't deviate too much.

This scheme is similar to that used by [NEMO-PISCES](https://www.nemo-ocean.eu/), although they scale the
tendency rather than the value, while other Earth system models simply set negative tracers to zero, for
example [NCAR's MARBL](https://marbl-ecosys.github.io/versions/latest_release/index.html) and
[NEMO-TOPAZ2](https://zenodo.org/records/2648099), which does not conserve mass. More complicated schemes
exist, for example [ROMS-BECS](https://zenodo.org/records/3988618) uses an implicite-itterative
approach where each component is updated in sequence to garantee mass conservation, possibly at the
expense of numerical precision.
"""
function ScaleNegativeTracers(tracers; scalefactors = ones(length(tracers)), invalid_fill_value = NaN, warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    tracers = Tuple(tracers)
    scalefactors = Tuple(scalefactors)
    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

"""
    ScaleNegativeTracers(bgc::AbstractBiogeochemistry; invalid_fill_value = NaN, warn = false)

Construct a negative-tracer treatment to scale the conserved tracers in `bgc` biogeochemistry.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(bgc::AbstractBiogeochemistry; invalid_fill_value = NaN, warn = false)
    tracers = conserved_tracers(bgc)

    return ScaleNegativeTracers(tracers; invalid_fill_value, warn)
end

function ScaleNegativeTracers(tracers::NTuple{<:Any, Symbol};
                              invalid_fill_value=NaN,
                              warn=false,)
    scalefactors = ntuple(_ -> 1.0, length(tracers))

    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

# multiple conserved groups
ScaleNegativeTracers(tracers::NamedTuple; invalid_fill_value = NaN, warn = false) =
    maybe_a_tuple(map(tn -> ScaleNegativeTracers(tn; invalid_fill_value, warn), values(tracers))...)

maybe_a_tuple(a) = a
maybe_a_tuple(args...) = tuple(args...)

function ScaleNegativeTracers(tracers::NamedTuple{<:Any, <:NTuple{<:Any, <:Number}};
                              invalid_fill_value=NaN,
                              warn=false,)
    scalefactors = values(tracers)
    tracer_names = keys(tracers)

    return ScaleNegativeTracers(tracer_names, scalefactors, invalid_fill_value, warn)
end

summary(scaler::ScaleNegativeTracers) = string("Mass conserving negative scaling of $(scaler.tracers)")
show(io::IO, scaler::ScaleNegativeTracers) = print(io, string(summary(scaler), "\n",
                                                          "└── Scalefactors: $(scaler.scalefactors)"))


function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)

    dev = device(architecture(model))

    apply_scale_for_negs!(dev, model, scale)

    return nothing
end

function apply_scale_for_negs!(dev::KA.GPU, model, scale)
    workgroup, worksize = work_layout(dev, model.grid, :xyz, (Center, Center, Center))

    scale_for_negs_kernel! = scale_for_negs_gpu!(dev, workgroup, worksize)

    field_scale = Tuple(Iterators.flatten(
        zip([model.tracers[tracer_name] for tracer_name in scale.tracers],
            scale.scalefactors
        )))

    scale_for_negs_kernel!(model.grid, scale.invalid_fill_value, field_scale...)

    return nothing
end

function apply_scale_for_negs!(dev::KA.CPU, model, scale)
    workgroup, worksize = work_layout(dev, model.grid, :xyz, (Center, Center, Center))

    scale_for_negs_kernel! = scale_for_negs_cpu!(dev, workgroup, worksize)

    # Fields themselves may have different types (due to type parameters)
    # It can cause type instability in the kernel if there is enough of different
    # field types in the tuple (and cause allocation in the for loop)
    # We need to homogenise types: hence we provide the data (aka OffsetArray)
    # directly
    tracers_to_scale = Tuple(model.tracers[tracer_name].data for tracer_name in scale.tracers)

    scale_for_negs_kernel!(model.grid, scale.invalid_fill_value, scale.scalefactors, tracers_to_scale)

    return nothing
end

# `CUDA.jl` has a tendency to make a local copy if we pass a list of tracer fields
# as a tuple (despite it being immutable). See the related issue:
# https://github.com/JuliaGPU/CUDA.jl/issues/1168
#
# This makes each thread use a large amount of 'Local' memory, which in turn
# introduces large amount of non-optimal memory traffic and slows the kernel
# significantly.
#
# The easiest workaround to avoid the local copy is to pass each field (and
# associated scalefactor) as a parameter to the kernel. But since we want to
# support arbitrary number of fields we need to make the kernel variadic.
#
# We expect Julia to inline the recursive calls and, effectively unroll the loops
@kernel cpu = false function scale_for_negs_gpu!(grid, invalid_fill_value, field_scale...)
    ijk = @index(Global, NTuple)

    if !immersed_cell(ijk..., grid)
        t, p = calculate_total_and_positive_part(0.0, 0.0, ijk, field_scale...)
    
        t = ifelse(t < 0, invalid_fill_value, t)
    
        correct_negative_fields!(t, p, ijk, field_scale...)
    end
    nothing
end

# Recursive step
@inline function calculate_total_and_positive_part(t,
                                                   p,
                                                   ijk,
                                                   field,
                                                   scale,
                                                   field_scale...)
    t, p = calculate_total_and_positive_part(t, p, ijk, field, scale)
    return calculate_total_and_positive_part(t, p, ijk, field_scale...)
end

# Recursion terminal
@inline function calculate_total_and_positive_part(t, p, ijk, field, scale)
    i, j, k = ijk
    value = @inbounds field[i, j, k]
    t += value * scale
    if value > 0
        p += value * scale
    end
    return t, p
end

# Recursive step
@inline function correct_negative_fields!(t, p, ijk, field, scale, field_scale...)
    correct_negative_fields!(t, p, ijk, field, scale)
    correct_negative_fields!(t, p, ijk, field_scale...)
    return nothing
end

# Recursion terminal
@inline function correct_negative_fields!(t, p, ijk, field, scale)
    i, j, k = ijk
    value = @inbounds field[i, j, k]
    new_value = ifelse(!isfinite(value) | (value > 0), value * t / p, 0)
    @inbounds field[i, j, k] = new_value
    return nothing
end

#
# The GPU kernel requires recursive implementation to avoid thread-local copies.
# However, when used on CPU it produces significant number of temporary allocations.
# This is most likely related to the fact that julia does not do tail call elimination
# on a CPU.
#
# Hence, for CPU code we need to fall-back to the loop-based version
#
@kernel function scale_for_negs_cpu!(grid, invalid_fill_value, scalefactors, fields)
    i, j, k = @index(Global, NTuple)

    if !immersed_cell(i, j, k, grid)
        t, p = 0.0, 0.0
    
        for (idx, field) in enumerate(fields)
            value = @inbounds field[i, j, k]
            scalefactor = @inbounds scalefactors[idx]
    
            t += value * scalefactor
            if value > 0
                p += value * scalefactor
            end
        end
    
        t = ifelse(t < 0, invalid_fill_value, t)
    
        for field in fields
            value = @inbounds field[i, j, k]
    
            new_value = ifelse(!isfinite(value) | (value > 0), value * t / p, 0)
    
            @inbounds field[i, j, k] = new_value
        end
    end
end
