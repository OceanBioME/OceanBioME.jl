# Individuals

The effects of individuals can be modelled in OceanBioME. We have implemented this through custom dynamics in the [Lagrangian Particle tracking feature of Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/model_setup/lagrangian_particles/). We have provided the tools to simulate spatially infinitesimal particles which can take nutrients etc. from a BGC model, and return waste. We also hope to provide models such as giant kelp which can act as a template for more complex individual behaviour.

To setup a particles first setup a particle struct as in Oceananigans with the desired properties, e.g.:
```@meta
DocTestSetup = quote
    using StructArrays, Roots
    using OceanBioME: Particles
end
```
``` jldoctest particles
struct Particle
    x :: AbstractFloat
    y :: AbstractFloat
    z :: AbstractFloat

    A :: AbstractFloat #amount particle has grown
    jₙ :: AbstractFloat #N uptake
    N :: AbstractFloat #nutrients availabel to the particle
end
n=2
particles = StructArray{Particle}((5*ones(n), 5*ones(n), zeros(n), ones(n), zeros(n), zeros(n)))
# output
2-element StructArray(::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}) with eltype Particle:
 Particle(5.0, 5.0, 0.0, 1.0, 0.0, 0.0)
 Particle(5.0, 5.0, 0.0, 1.0, 0.0, 0.0)

```

Then define how the particle evolves, in this fictitious case we will have a Mondo-quota nutrient uptake and growth:
``` jldoctest particles
function growth(A, N, params)
    jₙ = N/(params.K + N)
    
    return (A = jₙ, jₙ = jₙ)
end
# output
growth (generic function with 1 method)
```

We might also want to specify that instead of just being advected the particles do a random walk on the surface:
``` jldoctest particles
function custom_dynamics(particles, params)
    particles.properties.x .+= params.Δx*rand()
    particles.properties.y .+= params.Δy*rand()
    particles.properties.z .= 0.0
end
# output
custom_dynamics (generic function with 1 method)
```

We then use the provided setup function to turn this into particles for Oceananignans:
``` jldoctest particles
particles = Particles.setup(particles, growth, 
    (:A, :N), #arguments for growth
    (K = 1.0, Δx = 5.0, Δy = 5.0), #parameters
    (:A, ), #particles that should be integrated from the growth function output, i.e. +Δt * output
    (:jₙ, ), #particle properties updated by the growth function 
    (N = :N, ), #tracers being tracked, and the corresponding particle property to store in
    (N = (property=:jₙ, fallback=:N, fallback_scalefactor=1), ),#tracers being taken from, and the corresponding particle property saying how much to take
    10.0, #scale factor between particles and tracers (e.g. here 10*jₙ will be taken from N)
    custom_dynamics)
# output
2 LagrangianParticles with eltype Particle:
├── 6 properties: (:x, :y, :z, :A, :jₙ, :N)
├── particle-wall restitution coefficient: 1.0
├── 0 tracked fields: ()
└── dynamics: dynamics!

```
> WARNING: this API will imminently change to be more straightforward, hence the poor description of the sink fields

Please see the coming pages for specific particle models we have implimented (e.g. kelp)