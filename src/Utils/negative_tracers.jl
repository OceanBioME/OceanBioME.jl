using Oceananigans: fields, Simulation
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

import Adapt: adapt_structure

"""
    zero_negative_tracers!(sim; params = (exclude=(), warn=false))

Sets any tracers in `sim.model` which are negative to zero. Use like:
```julia
simulation.callbacks[:neg] = Callback(zero_negative_tracers!)
```
This is *NOT* a reccomended method to preserve positivity as it strongly does not conserve tracers.

Tracers to exclude can be set in the parameters and if `params.warn` is set to true a warning will be displayed when negative values are modified.
"""
function zero_negative_tracers!(sim; params = (exclude=(), warn=false))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if params.warn&&any(tracer .< 0.0) @warn "$tracer_name < 0" end
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

"""
    error_on_neg(sim; params = (exclude=(), ))

Throws an error if any tracers in `sim.model` are negative. Use like:
```julia
simulation.callbacks[:neg] = Callback(error_on_neg!)
```

Tracers to exclude can be set in the parameters.
"""
function error_on_neg!(sim; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if any(tracer .< 0.0) error("$tracer_name < 0") end
        end
    end
end

"""
    warn_on_neg(sim; params = (exclude=(), ))

Raises a warning if any tracers in `sim.model` are negative. Use like:
```julia
simulation.callbacks[:neg] = Callback(warn_on_neg!)
```

Tracers to exclude can be set in the parameters.
"""
function warn_on_neg!(sim; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if any(tracer .< 0.0) @warn "$tracer_name < 0" end
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

@kernel function scale_for_negs!(fields, tracers, scalefactors)
    i, j, k = @index(Global, NTuple)
    t, p = 0.0, 0.0
    @unroll for (i, tracer) in enumerate(tracers)
        field = fields[tracer]

        t += @inbounds field[i, j, k] * scalefactors[i]
        if field[i, j, k] > 0
            p += @inbounds field[i, j, k] * scalefactors[i]
        end
    end 
    @unroll for tracer in tracers
        field = fields[tracer]
        
        if @inbounds field[i, j, k] > 0
            @inbounds field[i, j, k] *= t / p
        else
            @inbounds field[i, j, k] = 0
        end
    end
end

struct ScaleNegativeTracers{FA, SA, W}
    tracers :: FA
    scalefactors :: SA
    warn :: W

    function ScaleNegativeTracers(tracers::FA, scalefactors::SA, warn::W) where {FA, SA, W}
        return new{FA, SA, W}(tracers, scalefactors, warn)
    end
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(adapt_structure(to, snt.tracers),
                                                                      adapt_structure(to, snt.scalefactors),
                                                                      adapt_structure(to, snt.warn))

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false)

Returns a callback that scales `tracers` so that none are negative. Use like:
```julia
negativity_protection! = ScaleNegativeTracers(tracers = (:P, :Z, :N))
simulation.callbacks[:neg] = Callback(negativity_protection!; callsite = UpdateStateCallsite())
```
This is a better but imperfect way to prevent numerical errors causing negative tracers. Please see discussion [here](https://github.com/OceanBioME/OceanBioME.jl/discussions/48). 
We plan to impliment positivity preserving timestepping in the future as the perfect alternative.
"""

function ScaleNegativeTracers(; model, tracers, scalefactors = NamedTuple{tracers}(ones(length(tracers))), warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    tracers = ntuple(n -> indexin([tracers[n]], [keys(model.tracers)...])[1], length(tracers))
    scalefactors = values(scalefactors)

    return ScaleNegativeTracers(tracers, scalefactors, warn)
end

@inline function (scale::ScaleNegativeTracers)(model)
    workgroup, worksize = work_layout(model.grid, :xyz)
    scale_for_negs_kernel! = scale_for_negs!(device(model.grid.architecture), workgroup, worksize)
    event = scale_for_negs_kernel!(model.tracers, scale.tracers, scale.scalefactors)
    wait(event)
end
@inline (scale::ScaleNegativeTracers)(sim::Simulation) = scale(sim.model) 

