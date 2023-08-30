using Oceananigans: fields, Simulation
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

import Adapt: adapt_structure, adapt
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!

"""
    ZeroNegativeTracers(; exclude = ())
Construct a modifier that zeros any negative tracers excluding those listed in `exclude`.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
@kwdef struct ZeroNegativeTracers{E}
    exclude :: E = ()
end

function update_biogeochemical_state!(model, zero::ZeroNegativeTracers)
    @unroll for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in zero.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{FA, SA, W}
    tracers :: FA
    scalefactors :: SA
    warn :: W

    ScaleNegativeTracers(tracers::FA, scalefactors::SA, warn::W) where {FA, SA, W} =
        new{FA, SA, W}(tracers, scalefactors, warn)
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(adapt(to, snt.tracers),
                                                                      adapt(to, snt.scalefactors),
                                                                      adapt(to, snt.warn))

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false)

Constructs a modifier to scale `tracers` so that none are negative. Use like:
```julia
modifier = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifier)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`zero_negative_tracers!`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

Future plans include implement a positivity-preserving timestepping scheme as the ideal alternative.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(tracers; scalefactors = NamedTuple{tracers}(ones(length(tracers))), warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    return ScaleNegativeTracers(tracers, scalefactors, warn)
end

"""
    ScaleNegativeTracers(model::UnderlyingBiogeochemicalModel; warn = false)

Constructs a modifier to scale the conserved tracers in `model`.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(model::UnderlyingBiogeochemicalModel; warn = false)
    tracers = conserved_tracers(model)
    scalefactors = NamedTuple{tracers}(ones(length(tracers)))

    return ScaleNegativeTracers(tracers, scalefactors, warn)
end

function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)
    workgroup, worksize = work_layout(model.grid, :xyz)
    dev = device(model.grid.architecture)
    scale_for_negs_kernel! = scale_for_negs!(dev, workgroup, worksize)
    scale_for_negs_kernel!(model.tracers, scale.tracers, scale.scalefactors)
end

@kernel function scale_for_negs!(fields, tracers, scalefactors)
    i, j, k = @index(Global, NTuple)
    t, p = 0.0, 0.0
    @unroll for (idx, tracer) in enumerate(tracers)
        field = @inbounds fields[tracer][i, j, k]
        scalefactor = @inbounds scalefactors[tracer]

        t += field * scalefactor
        if field > 0
            p += field * scalefactor
        end
    end 
    t < 0 && error("Cell total < 0, can not scale negative tracers.")
    @unroll for tracer in tracers
        field = @inbounds fields[tracer][i, j, k]
        
        if field > 0
            field *= t / p
        else
            field = 0
        end

        @inbounds fields[tracer][i, j, k] = field
    end
end