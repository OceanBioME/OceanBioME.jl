using Oceananigans: fields, Simulation
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device, architecture, on_architecture
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Adapt: adapt_structure, adapt
import Base: summary, show
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!

"""
    ZeroNegativeTracers(; exclude = ())

Construct a modifier that zeroes any negative tracers excluding those listed in `exclude`.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
@kwdef struct ZeroNegativeTracers{E}
    exclude :: E = ()
end

function update_biogeochemical_state!(model, zero::ZeroNegativeTracers)
    for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in zero.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{FA, SA, FV, W}
           tracers :: FA
      scalefactors :: SA
invalid_fill_value :: FV
              warn :: W

    ScaleNegativeTracers(tracers::FA, scalefactors::SA, invalid_fill_value::FV, warn::W) where {FA, SA, W, FV} =
        warn ? error("Warning not currently implemented") : new{FA, SA, FV, W}(tracers, scalefactors, invalid_fill_value, warn)
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(adapt(to, snt.tracers),
                                                                      adapt(to, snt.scalefactors),
                                                                      adapt(to, snt.invalid_fill_value),
                                                                      adapt(to, snt.warn))

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false, invalid_fill_value = NaN)

Constructs a modifier to scale `tracers` so that none are negative. Use like:
```julia
modifier = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifier)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`ZeroNegativeTracers`](@ref). Please see [discussion in
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

    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

"""
    ScaleNegativeTracers(bgc::AbstractBiogeochemistry; warn = false)

Construct a modifier to scale the conserved tracers in `bgc` biogeochemistry.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(bgc::AbstractBiogeochemistry, grid; invalid_fill_value = NaN, warn = false)
    tracers = conserved_tracers(bgc)

    return ScaleNegativeTracers(tracers, grid; invalid_fill_value, warn)
end

# for when `conserved_tracers` just returns a tuple of symbols
function ScaleNegativeTracers(tracers, grid; invalid_fill_value=NaN, warn=false)
    scalefactors = ones(length(tracers))
    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

function ScaleNegativeTracers(tracers::NTuple{<:Any,Symbol},
                              grid;
                              invalid_fill_value=NaN,
                              warn=false,)
    scalefactors = ones(length(tracers))

    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

# multiple conserved groups
ScaleNegativeTracers(tracers::Tuple, grid;invalid_fill_value = NaN, warn = false) =
    tuple(map(tn -> ScaleNegativeTracers(tn, grid; invalid_fill_value, warn), tracers)...)

function ScaleNegativeTracers(tracers::NamedTuple,
                              grid;
                              invalid_fill_value=NaN,
                              warn=false,)
    scalefactors = [tracers.scalefactors...]
    tracer_names = tracers.tracers

    return ScaleNegativeTracers(tracer_names, scalefactors, invalid_fill_value, warn)
end

summary(scaler::ScaleNegativeTracers) = string("Mass conserving negative scaling of $(scaler.tracers)")
show(io::IO, scaler::ScaleNegativeTracers) = print(io, string(summary(scaler), "\n",
                                                          "└── Scalefactors: $(scaler.scalefactors)"))

# Helper function to create interleaved tuple from two lists
# is not very robust and supports only lists of equal length!
interleave(l_list, r_list) = interleave_step(r_list, l_list...)
function interleave_step(l_list, r_start, r_tail...)
    return (r_start, interleave_step(r_tail, l_list...)...)
end
interleave_step(r_start) = r_start

function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)
    workgroup, worksize = work_layout(model.grid, :xyz, (Center, Center, Center))

    dev = device(architecture(model))

    scale_for_negs_kernel! = scale_for_negs!(dev, workgroup, worksize)

    tracers_to_scale = Tuple(model.tracers[tracer_name] for tracer_name in scale.tracers)

    field_scale = interleave([model.tracers[tracer_name] for tracer_name in scale.tracers],
                             scale.scalefactors)

    scale_for_negs_kernel!(scale.invalid_fill_value, field_scale...)

    return nothing
end

# `CUDA.jl` has a tendency to make a local copy if we pass a list of tracer fields
# as a tuple (despite it beeing immutable). See the related issue:
# https://github.com/JuliaGPU/CUDA.jl/issues/1168
#
# This makes each thread use a large amout of 'Local' memory, which in turn 
# introduces large amout of unoptimal memory traffic and slows the kernel
# singnificantly.
#
# The easiest workaround to avoid the local copy is to pass each field (and
# associated scalfactor) as a parameter to the kernel. But since we want to
# support arbitrary number of fields we need to make the kernel variadic.
#
# We expect Julia to inline the recursive calls and, effectivly unroll the loops
@kernel function scale_for_negs!(invalid_fill_value, field_scale...)
    ijk = @index(Global, NTuple)

    t, p = reduce_over_fields(0.0, 0.0, ijk, field_scale...)

    t = ifelse(t < 0, invalid_fill_value, t)

    update_fields!(t, p, ijk, field_scale...)
    nothing
end

# Recursive step
@inline function reduce_over_fields(t,
                                    p,
                                    ijk,
                                    field::T,
                                    scale::Real,
                                    field_scale...) where {T}
    t, p = reduce_over_fields(t, p, ijk, field, scale)
    return reduce_over_fields(t, p, ijk, field_scale...)
end

# Recursion terminal
@inline function reduce_over_fields(t, p, ijk, field::T, scale::Real) where {T}
    i, j, k = ijk
    value = @inbounds field[i, j, k]
    t += value * scale
    if value > 0
        p += value * scale
    end
    return t, p
end

# Recursive step
@inline function update_fields!(t, p, ijk, field::T, scale::Real, field_scale...) where {T}
    update_fields!(t, p, ijk, field, scale)
    return update_fields!(t, p, ijk, field_scale...)
end

# Recursion terminal
@inline function update_fields!(t, p, ijk, field::T, scale::Real) where {T}
    i, j, k = ijk
    value = @inbounds field[i, j, k]
    new_value = ifelse(!isfinite(value) | (value > 0), value * t / p, 0)
    @inbounds field[i, j, k] = new_value
end
