using Oceananigans.Architectures: on_architecture, device, architecture
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_auxiliary_fields
using Oceananigans.Models.LagrangianParticleTracking: ParticleBoundaryConditions, _advect_particles!, total_velocities

# put nothing to do nothing
advect_particles!(advection, particles, model, Δt) = nothing

"""
    LagrangianAdvection

Specifies that particles should move in a purley lagrangian mannor.
"""
struct LagrangianAdvection end

function advect_particles!(::LagrangianAdvection, particles, model, Δt)
    workgroup = min(length(particles), 256)
    worksize = length(particles)

    arch = architecture(model)

    # Advect particles
    advect_particles_kernel! = _advect_particles!(device(arch), workgroup, worksize)

    # Reflecting boundaries with a restitution of 1, and no absorption record since no particle is removed
    advect_particles_kernel!(particles, ParticleBoundaryConditions(), nothing, model.grid, Δt, total_velocities(model))

    return nothing
end
