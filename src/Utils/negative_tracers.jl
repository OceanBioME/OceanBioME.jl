using Oceananigans: fields, Simulation
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

import Adapt: adapt_structure, adapt

"""
    zero_negative_tracers!(sim; params = (exclude=(), ))

Sets any tracers in `sim.model` which are negative to zero. Use like:
```julia
simulation.callbacks[:neg] = Callback(zero_negative_tracers!)
```

Tracers to exclude can be set in the parameters.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
function zero_negative_tracers!(model; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in params.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end
@inline zero_negative_tracers!(sim::Simulation; params = (exclude=(), )) = zero_negative_tracers!(sim.model; params)

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
    @unroll for (idx, tracer) in enumerate(tracers)
        field = @inbounds fields[tracer][i, j, k]
        scalefactor = @inbounds scalefactors[idx]

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

Returns a callback that scales `tracers` so that none are negative. Use like:
```julia
negativity_protection! = ScaleNegativeTracers(; model, tracers = (:P, :Z, :N))
simulation.callbacks[:neg] = Callback(negativity_protection!; callsite = UpdateStateCallsite())
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`zero_negative_tracers!`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

Future plansj include implement a positivity=preserving timestepping scheme as the ideal alternative.
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
    dev = device(model.grid.architecture)
    scale_for_negs_kernel! = scale_for_negs!(dev, workgroup, worksize)
    scale_for_negs_kernel!(values(model.tracers), scale.tracers, scale.scalefactors)
end
@inline (scale::ScaleNegativeTracers)(sim::Simulation) = scale(sim.model) 

@kernel function _remove_NaN_tendencies!(fields)
    i, j, k = @index(Global, NTuple)
    for field in fields
        if @inbounds isnan(field[i, j, k])
            @inbounds field[i, j, k] = 0.0
        end
    end
end

"""
    remove_NaN_tendencies!(model)

Zeros any `NaN` value tendencies as a final protection against negative tracer run away.
"""
@inline function remove_NaN_tendencies!(model)
    workgroup, worksize = work_layout(model.grid, :xyz)
    remove_NaN_tendencies_kernel! = _remove_NaN_tendencies!(device(model.grid.architecture), workgroup, worksize)
    remove_NaN_tendencies_kernel!(values(model.timestepper.Gⁿ))

    return nothing
end
