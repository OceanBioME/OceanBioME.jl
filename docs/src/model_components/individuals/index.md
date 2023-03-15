# [Individuals](@id individuals)

The effects of individuals can be modelled in OceanBioME. We have implemented this through custom dynamics in the [Lagrangian Particle tracking feature of Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/model_setup/lagrangian_particles/). We have extended these functionalities to make it easier to implement "active" particles which interact with the tracers. We have then implemented a model of [sugar kelp](@ref SLatissima) which can be followed as an example of using this functionality.

To setup particles first create a particle struct with the desired properties, e.g.:
```@meta
DocTestSetup = quote
    using OceanBioME.Particles: BiogeochemicalParticles, get_node
    using Oceananigans.Fields: interpolate
end
```
``` jldoctest particles
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

``` jldoctest particles
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.LagrangianParticleTracking: update_particle_properties!, _advect_particles!
```

First, to integrate the particles properties we overload `update_particle_properties`, in this fictitious case we will have a Mondo-quota nutrient uptake and growth:

``` jldoctest particles

function update_particle_properties!(particles::GrowingParticles, model, bgc, Δt)
    @inbounds for p in 1:length(particles)
        nutrients = @inbounds interpolate(model.tracers.NO₃, particle.x[p], particle.y[p], particle.z[p])

        uptake = nutrients / (particle.nutrients_half_saturation + nutrients)

        particles.size[p] += uptake * Δt
        particles.nitrate_uptake[p] = uptake
    end
end
```

In this example the particles will not move around, and are only integrated on a single thread. For a more comprehensive example see the [Sugar Kelp](@ref SLatissima) implementation. We then need to update the tracer tendencies to match the nutrients' uptake:

``` jldoctest particles

function update_tendencies!(bgc, particles::GrowingParticles, model)
    @inbounds for p in 1:length(particles)
        # here we use an OceanBioME utility to find the nearest nodes to apply the tendency to
        nodes, normfactor = @inbounds get_nearest_nodes(particles.x[p], particles.y[p], particles.z[p], model.grid, (Center(), Center(), Center()))

        for (i, j, k, d) in nodes 
            # Reflect back on Bounded boundaries or wrap around for Periodic boundaries
            i, j, k = (get_node(TX(), i, grid.Nx), get_node(TY(), j, grid.Ny), get_node(TZ(), k, grid.Nz))

            node_volume = volume(i, j, k, grid, LX(), LY(), LZ())
            @inbounds model.timestepper.Gⁿ.NO₃[i, j, k] += particles.nitrate_uptake[p] / (d * node_volume)
        end
    end
end
```

Now we can just plug this into any biogeochemical model setup to have particles (currently [NPZD](@ref NPZD) and [LOBSTER](@ref LOBSTER)):
``` jldoctest particles
# Start the particles randomly distributed, floating on the surface
Lx, Ly, Lz = 1000, 1000, 100

grid = RectilinearGrid(; size = (64, 64, 16), extent = (Lx, Ly, Lz))

particles = GrowingParticles(0.5, zeros(3), zeros(3), rand(3) * Lx, rand(3) * Ly, zeros(3))

biogeochemistry = LOBSTER(; grid, particles)
```