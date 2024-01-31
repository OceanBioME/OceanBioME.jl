# [Individuals](@id individuals)

The effects of individuals can be modelled in OceanBioME. We have implemented this through custom dynamics in the [Lagrangian Particle tracking feature of Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/model_setup/lagrangian_particles/). We have extended these functionalities to make it easier to implement "active" particles which interact with the tracers. We have then implemented a model of [sugar kelp](@ref SLatissima) which can be followed as an example of using this functionality.

To setup particles first create a particle type with the desired properties, e.g.:

```@example particles
using OceanBioME.Particles: BiogeochemicalParticles

struct GrowingParticles{FT, VT} <: BiogeochemicalParticles 
    nutrients_half_saturation :: FT

    size :: VT
    nitrate_uptake :: VT

    x :: VT
    y :: VT
    z :: VT
end
```

You then need to overload particular functions to integrate the growth, so they need to first be `import`ed:

```@example particles
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Models.LagrangianParticleTracking: update_lagrangian_particle_properties!
```

First, to integrate the particles properties we overload `update_lagrangian_particle_properties!`;
in this fictitious case we will have a Mondo-quota nutrient uptake and growth:

```@example particles
using Oceananigans.Fields: interpolate

function update_lagrangian_particle_properties!(particles::GrowingParticles, model, bgc, Δt)
    @inbounds for p in 1:length(particles)
        nutrients = @inbounds interpolate(model.tracers.NO₃, particle.x[p], particle.y[p], particle.z[p])

        uptake = nutrients / (particle.nutrients_half_saturation + nutrients)

        particles.size[p] += uptake * Δt
        particles.nitrate_uptake[p] = uptake
    end
    return nothing
end

nothing #hide
```

In this example the particles will not move around, and are only integrated on a single thread. For a more comprehensive example see the [Sugar Kelp](@ref SLatissima) implementation. We then need to update the tracer tendencies to match the nutrients' uptake:

```@example particles
using OceanBioME.Particles: get_node

function update_tendencies!(bgc, particles::GrowingParticles, model)
    @inbounds for p in 1:length(particles)
        i, j, k = fractional_indices((x, y, z), grid, Center(), Center(), Center())

        # Convert fractional indices to unit cell coordinates 0 ≤ (ξ, η, ζ) ≤ 1
        # and integer indices (with 0-based indexing).
        ξ, i = modf(i)
        η, j = modf(j)
        ζ, k = modf(k)

        # Round to nearest node and enforce boundary conditions
        i, j, k = (get_node(TX(), Int(ifelse(ξ < 0.5, i + 1, i + 2)), grid.Nx),
                   get_node(TY(), Int(ifelse(η < 0.5, j + 1, j + 2)), grid.Ny),
                   get_node(TZ(), Int(ifelse(ζ < 0.5, k + 1, k + 2)), grid.Nz))

        node_volume = volume(i, j, k, grid, Center(), Center(), Center())

        model.timestepper.Gⁿ.NO₃[i, j, k] += particles.nitrate_uptake[p] / (d * node_volume)
    end
    return nothing
end

nothing #hide
```

Now we can just plug this into any biogeochemical model setup to have particles (currently [NPZD](@ref NPZD) and [LOBSTER](@ref LOBSTER)):

```@example particles
using OceanBioME, Oceananigans

Lx, Ly, Lz = 1000, 1000, 100
grid = RectilinearGrid(; size = (64, 64, 16), extent = (Lx, Ly, Lz))

# Start the particles randomly distributed, floating on the surface
particles = GrowingParticles(0.5, zeros(3), zeros(3), rand(3) * Lx, rand(3) * Ly, zeros(3))

biogeochemistry = LOBSTER(; grid, particles)
```
