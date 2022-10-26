using Oceananigans: fields
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device

function zero_negative_tracers!(sim; params = (exclude=(), warn=false))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if params.warn&&any(tracer .< 0.0) @warn "$tracer_name < 0" end
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

function error_on_neg!(sim; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if any(tracer .< 0.0) error("$tracer_name < 0") end
        end
    end
end

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

function scale_negative_tracers!(sim, params=(conserved_group = (), warn=false)) #this can be used to conserve sub groups e.g. just saying NO₃ and NH₄ 
    workgroup, worksize = work_layout(sim.model.grid, :xyz)
    scale_for_negs_kernel! = scale_for_negs!(device(sim.model.grid.architecture), workgroup, worksize)
    model_fields = fields(sim.model)
    event = scale_for_negs_kernel!(model_fields[params.conserved_group], params.warn)
    wait(event)
end