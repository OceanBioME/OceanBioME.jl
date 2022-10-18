using KernelAbstractions.Extras.LoopInfo: @unroll

function no_negative_tracers!(sim; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            if any(tracer .< 0.0) @warn "$tracer_name < 0" end
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end