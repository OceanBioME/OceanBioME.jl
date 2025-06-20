module DroguedBuoyModel

export DroguedBuoyDynamics

using Oceananigans.Architectures: architecture
using Oceananigans.Models: total_velocities
using Oceananigans.Models.LagrangianParticleTracking: LagrangianParticles, advect_particle
using Oceananigans.Utils: KernelParameters, launch!
using KernelAbstractions: @index, @kernel

import Oceananigans.Models.LagrangianParticleTracking: advect_lagrangian_particles!

struct DroguedBuoyDynamics{DD}
    drogue_depths :: DD
end

@inline (::DroguedBuoyDynamics)(args...) = nothing

const DroguedBuoyParticle = LagrangianParticles{<:Any, <:Any, <:Any, <:DroguedBuoyDynamics}

function advect_lagrangian_particles!(particles::DroguedBuoyParticle, model, Δt)
    grid = model.grid
    arch = architecture(grid)
    parameters = KernelParameters(1:length(particles))

    launch!(arch, grid, parameters,
            _advect_drogued_particles!,
            particles.properties, particles.restitution, model.grid, Δt, total_velocities(model), particles.dynamics.drogue_depths)

    return nothing
end

@kernel function _advect_drogued_particles!(properties, restitution, grid, Δt, velocities, drogue_depths)
    p = @index(Global)

    @inbounds begin
        x = properties.x[p]
        y = properties.y[p]
        z = drogue_depths[p]
    end

    x⁺, y⁺, _ = advect_particle((x, y, z), properties, p, restitution, grid, Δt, velocities)

    @inbounds begin
        properties.x[p] = x⁺
        properties.y[p] = y⁺
    end
end

end # module