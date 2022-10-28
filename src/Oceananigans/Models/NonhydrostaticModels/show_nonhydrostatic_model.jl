using Oceananigans.Utils: prettytime, ordered_dict_show, prettykeys
using Oceananigans.TurbulenceClosures: closure_summary

import Base: show

function show(io::IO, model::NonhydrostaticModel)
    TS = nameof(typeof(model.timestepper))
    tracernames = prettykeys(model.tracers)
    
    print(io, summary(model), "\n",
        "├── grid: ", summary(model.grid), "\n",
        "├── timestepper: ", TS, "\n",
        "├── tracers: ", tracernames, "\n",
        "├── closure: ", closure_summary(model.closure), "\n",
        "├── buoyancy: ", summary(model.buoyancy), "\n")

    if !isnothing(model.auxiliary_fields)
        print(io, "├── auxiliary_fields: ", prettykeys(model.auxiliary_fields), "\n")
    end

    if isnothing(model.particles)
        print(io, "└── coriolis: ", summary(model.coriolis))
    else
        particles = model.particles.properties
        properties = propertynames(particles)
        print(io, "├── coriolis: ", summary(model.coriolis), "\n")
        print(io, "└── particles: ", summary(model.particles))
    end
end

