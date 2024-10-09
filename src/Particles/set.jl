using Oceananigans.Architectures: architecture, on_architecture

import Oceananigans.Fields: set!

function set!(particles::BiogeochemicalParticles; kwargs...)
    arch = architecture(particles)

    for (n, v) in pairs(kwargs)
        if n == :x
            particles.x .= on_architecture(arch, v)
        elseif n == :y
            particles.y .= on_architecture(arch, v)
        elseif n == :z
            particles.z .= on_architecture(arch, v)
        elseif n == :scalefactors
            particles.scalefactors .= on_architecture(arch, v)
        else
            particles.fields[n] .= on_architecture(arch, v)
        end
    end
    
    return nothing
end

