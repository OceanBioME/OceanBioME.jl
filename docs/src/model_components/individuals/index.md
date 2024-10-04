# [Individuals](@id individuals)

The effects of individuals can be modelled in OceanBioME. We have implemented this through custom dynamics in the [Lagrangian Particle tracking feature of Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/model_setup/lagrangian_particles/). We have extended these functionalities to make it easier to implement "active" particles which interact with the tracers. We have then implemented a model of [sugar kelp](@ref sugar-kelp) which can be followed as an example of using this functionality.

To setup particles first create a particle biogeochemistry, e.g.:

```@example particles

struct GrowingParticles{FT}
    nutrients_half_saturation :: FT
end
```

We then need to add some methods to tell `OceanBioME` what properties this particle has, and what tracers it interacts with:

```@example particles

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers

required_particle_fields(::GrowingParticles) = (:S, )
required_tracers(::GrowingParticles) = (:N, )
coupled_tracers(::GrowingParticles) = (:N, )

```

So our model is going to track the `S`ize of the particles and take up `N`utrients. 
Now we need to how this growth happens. 
The forcing functions should be of the form `(particles::ParticleBiogeochemistry)(::Val{:PROPERTY}, t, required_particle_fields..., required_tracers...)`, so in this example:
```@example particles
(p::GrowingParticles)(::Val{:S}, t, S, N) = N / (N + p.nutrient_half_saturation)
(p::GrowingParticles)(::Val{:N}, t, S, N) = - N / (N + p.nutrient_half_saturation)
```

We can then create an instance of this particle model using `BiogeochemicalParticles`, and set their initial position and size:
```@example particles
using OceanBioME, Oceananigans

Lx, Ly, Lz = 100, 100, 100
grid = RectilinearGrid(; size = (8, 8, 8), extent = (Lx, Ly, Lz))

particles = BiogeochemicalParticles(10; grid, biogeochemistry = GrowingParticles())

set!(particles, S = 0.1, x = rand(10) * Lx, y = rand(10) * Ly, z = rand(10) * Lz)
```

We can then put these into a compatible biogeochemical model, for example:
```@example particles

biogeochemistry = NPZD(; grid, particles)
```
