using KernelAbstractions.Extras.LoopInfo: @unroll

function no_negative_tracers!(sim; params = (exclude=(), ))
    @unroll for (tracer_name, tracer) in pairs(sim.model.tracers)
        if !(tracer_name in params.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end