using Oceananigans: fields, Simulation
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

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

@kernel function scale_for_negs!(fields, warn)
    i, j, k = @index(Global, NTuple)
    t, p = 0.0, 0.0
    @unroll for field in fields
        t += @inbounds field[i, j, k]
        if field[i, j, k] > 0
            p += @inbounds field[i, j, k]
        end
    end 
    @unroll for field in fields
        if @inbounds field[i, j, k]>0
            @inbounds field[i, j, k] *= t/p
        else
            if warn @warn "Scaling negative" end
            @inbounds field[i, j, k] = 0
        end
    end
end

"""
    scale_negative_tracers!(sim, params=(conserved_group = (), warn=false))

Scales tracers in `conserved_group` so that none are negative. Use like:
```julia
simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))
```
This is a better but imperfect way to prevent numerical errors causing negative tracers. Please see discussion [here](https://github.com/OceanBioME/OceanBioME.jl/discussions/48). 
We plan to impliment positivity preserving timestepping in the future as the perfect alternative.

Tracers conserve should be set in the parameters and if `params.warn` is set to true a warning will be displayed when negative values are modified.
"""
function scale_negative_tracers!(model, params=(conserved_group = (), warn=false)) #this can be used to conserve sub groups e.g. just saying NO₃ and NH₄ 
    workgroup, worksize = work_layout(model.grid, :xyz)
    scale_for_negs_kernel! = scale_for_negs!(device(model.grid.architecture), workgroup, worksize)
    model_fields = fields(model)
    event = scale_for_negs_kernel!(model_fields[params.conserved_group], params.warn)
    wait(event)
end
@inline scale_negative_tracers!(sim::Simulation, args...) = scale_negative_tracers!(sim.model, args...) 