"""
This module adds active particles which sink, dissolve, and release material into tracer fields. For example, mineral particles added to the ocean for the purpose of ocean alkalinity enhancement (OAE) or dissolving secondarary mineral precipitate.

The particles are defined by a radius, r (m), which sets how fast they dissolve and sink.

Prognostic properties
=====================
* Radius: r (m)

Tracer dependencies
===================
Those listed in `required_fields`, by default none

Coupled tracers
===============
Those listed in `released_tracers`, by default alkalinity: Alk (mmol/m³)
"""

module DissolvingParticlesModel

export DissolvingMineral, DissolvingParticles, StokesSinking

using Adapt

using Oceananigans: fields
using Oceananigans.Grids: AbstractGrid

using OceanBioME.Particles: BiogeochemicalParticles, LagrangianAdvection, extract_tracer_values

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers, update_particle_state!, advect_particles!
import Oceananigans.Models.LagrangianParticleTracking: particle_w_velocity

import Base: summary, show

"""
    StokesSinking(FT = Float64;
                  particle_density,
                  water_density = 1025,
                  dynamic_viscosity = 1.08e-3,
                  gravitational_acceleration = 9.81)

Function to calculate Stokes' law for the speed (m/s) at which a sphere of radius `r` (m) sinks,
called as `stokes(r, field_values...)`:

```math
w = \\frac{2(\\rho_p - \\rho)gr^2}{9\\mu},
```

where ``\\rho_p`` is the `particle_density` (kg/m³), ``\\rho`` the `water_density` (kg/m³),
``\\mu`` the `dynamic_viscosity` (kg/m/s), and ``g`` the `gravitational_acceleration` (m/s²).

The `water_density` can be a number or a function of the values of the mineral's `required_fields`
(see [`DissolvingMineral`](@ref)), e.g. `water_density = (T, S) -> ...` with `required_fields = (:T, :S)`.

This is only valid when the particle Reynolds number, ``2rw\\rho/\\mu``, is small,
i.e. for particles smaller than about 100 μm.
"""
struct StokesSinking{FT, R}
              particle_density :: FT
                 water_density :: R
             dynamic_viscosity :: FT
    gravitational_acceleration :: FT
end

function StokesSinking(FT = Float64;
                       particle_density,
                       water_density = 1025,
                       dynamic_viscosity = 1.08e-3,
                       gravitational_acceleration = 9.81)

    water_density = maybe_convert(FT, water_density)

    return StokesSinking{FT, typeof(water_density)}(particle_density, water_density, dynamic_viscosity, gravitational_acceleration)
end

@inline function (stokes::StokesSinking)(r, field_values...)
    ρₚ = stokes.particle_density
    ρ  = evaluate(stokes.water_density, field_values...)
    μ  = stokes.dynamic_viscosity
    g  = stokes.gravitational_acceleration

    return 2 * (ρₚ - ρ) * g * r^2 / (9 * μ)
end

"""
    DissolvingMineral(FT = Float64;
                      density = 2211,
                      molar_mass = 74.093e-3,
                      released_tracers = (; Alk = 2),
                      dissolution_rate = 1e-8,
                      sinking_speed = StokesSinking(FT; particle_density = density),
                      required_fields = ())

Mineral particles which dissolve and release material into tracer fields as they sink.
The only property tracked is the particle radius, `r` (m), which sets how fast the particles dissolve and sink.

The radius shrinks at the `dissolution_rate`,

```math
\\frac{dr}{dt} = -k(r),
```

so each particle loses volume at a rate of ``4\\pi r^2 k(r)`` and releases tracers at
``4\\pi r^2 k(r) \\nu / V_m``, where ``\\nu`` is the moles of tracer released per mole of mineral (`released_tracers`) and ``V_m`` is the molar volume (`molar_mass / density`).
The released tracers are converted from moles to mmol to match the units of the tracers in OceanBioME.

Particles also sink relative to the water at `sinking_speed` when they are advected by
`LagrangianAdvection` (the default for `BiogeochemicalParticles`).

`dissolution_rate` and `sinking_speed` can depend on model fields (tracers, velocities, or auxiliary fields)
at the particle position by listing their names in `required_fields`. The functions are then called with the
radius followed by the field values, e.g. `required_fields = (:T, )` and `dissolution_rate = (r, T) -> ...`.

Keyword Arguments
=================

- `density`: density of the mineral (kg/m³)
- `molar_mass`: molar mass of the mineral (kg/mol)
- `released_tracers`: `NamedTuple` of the moles of each tracer released per mole of mineral dissolved,
  e.g. `(; Alk = 2, DIC = 1)` for calcium carbonate
- `dissolution_rate`: the speed (m/s) that the radius shrinks, ``k(r)``, either a number or a function of radius (and `required_fields`)
- `sinking_speed`: the speed (m/s, positive downwards) that the particles sink, either a number or a function of radius (and `required_fields`)
- `required_fields`: names of the model fields that `dissolution_rate` and `sinking_speed` depend on

The default values of `density`, `molar_mass`, and `released_tracers` are for calcium hydroxide (hydrated lime, Ca(OH)₂), and by default the particles sink according to Stokes' law (see [`StokesSinking`](@ref)).
The default `dissolution_rate` of 1e-8 m/s is only illustrative (a 100 μm particle dissolves in about 3 hours), and should be set for the mineral and conditions being modelled.

Example
=======

```jldoctest
julia> using OceanBioME

julia> calcite = DissolvingMineral(density = 2710, molar_mass = 100.09e-3, released_tracers = (; Alk = 2, DIC = 1))
DissolvingMineral{Float64} biogeochemistry tracking particle `r`adius
├── density: 2710.0 kg/m³
├── molar mass: 0.10009 kg/mol
├── released tracers: (Alk = 2.0, DIC = 1.0) mol/mol
├── dissolution rate: 1.0e-8 m/s
├── sinking speed: StokesSinking{Float64, Float64}(2710.0, 1025.0, 0.00108, 9.81) m/s
└── required fields: ()
```
"""
struct DissolvingMineral{FT, R, D, W, F}
             density :: FT
          molar_mass :: FT
    released_tracers :: R
    dissolution_rate :: D
       sinking_speed :: W

    # the names of the `required_fields` are a type parameter so that they are known inside kernels
    DissolvingMineral(density::FT, molar_mass::FT, released_tracers::R, dissolution_rate::D, sinking_speed::W,
                      required_fields::Tuple) where {FT, R, D, W} =
        new{FT, R, D, W, required_fields}(density, molar_mass, released_tracers, dissolution_rate, sinking_speed)
end

function DissolvingMineral(FT = Float64;
                           density = 2211,
                           molar_mass = 74.093e-3,
                           released_tracers = (; Alk = 2),
                           dissolution_rate = 1e-8,
                           sinking_speed = StokesSinking(FT; particle_density = density),
                           required_fields = ())

    released_tracers = map(ν -> convert(FT, ν), released_tracers)
    dissolution_rate = maybe_convert(FT, dissolution_rate)
    sinking_speed = maybe_convert(FT, sinking_speed)

    return DissolvingMineral(convert(FT, density), convert(FT, molar_mass), released_tracers, dissolution_rate, sinking_speed,
                             required_fields)
end

maybe_convert(FT, value::Number) = convert(FT, value)
maybe_convert(FT, func) = func

# `dissolution_rate`, `sinking_speed`, and `water_density` can either be numbers or functions
@inline evaluate(value::Number, args...) = value
@inline evaluate(func, args...) = func(args...)

"""
    DissolvingParticles(n; grid, mineral_parameters = NamedTuple(), kwargs...)

Sets up `n` dissolving `BiogeochemicalParticles` with default [`DissolvingMineral`](@ref) parameters except those specified in `mineral_parameters`. `kwargs` are passed onto `BiogeochemicalParticles`.
"""
DissolvingParticles(n; grid::AbstractGrid{FT}, mineral_parameters = NamedTuple(), kwargs...) where FT =
    BiogeochemicalParticles(n; grid, biogeochemistry = DissolvingMineral(FT; mineral_parameters...), kwargs...)

@inline required_particle_fields(::DissolvingMineral) = (:r, )
@inline required_tracers(::DissolvingMineral{<:Any, <:Any, <:Any, <:Any, F}) where F = F
@inline coupled_tracers(mineral::DissolvingMineral) = keys(mineral.released_tracers)

# Volume of mineral dissolving per second (m³/s), which stops once the particle has fully dissolved
@inline dissolution(mineral, r, field_values...) =
    ifelse(r > 0, 4π * r^2 * evaluate(mineral.dissolution_rate, r, field_values...), zero(r))

# Change in radius (m/s)
@inline (mineral::DissolvingMineral)(::Val{:r}, t, r, field_values...) =
    - ifelse(r > 0, evaluate(mineral.dissolution_rate, r, field_values...), zero(r))

# Release of each tracer (mmol/s)
@inline function (mineral::DissolvingMineral)(::Val{tracer}, t, r, field_values...) where tracer
    molar_volume = mineral.molar_mass / mineral.density # m³/mol
    ν = mineral.released_tracers[tracer] # mol tracer / mol mineral

    return 1000 * ν * dissolution(mineral, r, field_values...) / molar_volume
end

@inline sinking_speed(mineral, r, field_values...) =
    ifelse(r > 0, evaluate(mineral.sinking_speed, r, field_values...), zero(r))

const DissolvingBiogeochemicalParticles = BiogeochemicalParticles{<:Any, <:DissolvingMineral}

# Oceananigans' particle advection only passes the particles into `particle_w_velocity`, so we
# wrap them with the model fields so that the sinking speed can depend on the `required_fields`
struct ParticlesWithFields{P, G, F, X, Y, Z}
    particles :: P
         grid :: G
       fields :: F
            x :: X # read by Oceananigans' particle advection
            y :: Y
            z :: Z
end

ParticlesWithFields(particles, grid, fields) =
    ParticlesWithFields(particles, grid, fields, particles.x, particles.y, particles.z)

Adapt.adapt_structure(to, p::ParticlesWithFields) =
    ParticlesWithFields(adapt(to, p.particles), adapt(to, p.grid), adapt(to, p.fields))

Base.length(p::ParticlesWithFields) = length(p.particles)

advect_particles!(advection::LagrangianAdvection, particles::DissolvingBiogeochemicalParticles, model, Δt) =
    advect_particles!(advection, ParticlesWithFields(particles, model.grid, fields(model)), model, Δt)

# The particles sink relative to the water when they're moved by `LagrangianAdvection`
@inline function particle_w_velocity(p::ParticlesWithFields, n, w)
    particles = p.particles

    r = @inbounds particles.fields.r[n]
    field_values = extract_tracer_values(nothing, particles.field_interpolation, particles, p.grid, p.fields, n)

    return w - sinking_speed(particles.biogeochemistry, r, field_values...)
end

# The final timestep of dissolution can overshoot to a negative radius, so we clip it
# (this doesn't change what was released because nothing is released when `r ≤ 0`)
function update_particle_state!(particles::DissolvingBiogeochemicalParticles, model, Δt)
    r = particles.fields.r
    r .= max.(r, 0)
    return nothing
end

summary(::DissolvingMineral{FT}) where FT = "DissolvingMineral{$FT} biogeochemistry"

show(io::IO, mineral::DissolvingMineral) =
    print(io, summary(mineral), " tracking particle `r`adius", "\n",
              "├── density: ", mineral.density, " kg/m³", "\n",
              "├── molar mass: ", mineral.molar_mass, " kg/mol", "\n",
              "├── released tracers: ", mineral.released_tracers, " mol/mol", "\n",
              "├── dissolution rate: ", mineral.dissolution_rate, " m/s", "\n",
              "├── sinking speed: ", mineral.sinking_speed, " m/s", "\n",
              "└── required fields: ", required_tracers(mineral))

end # module
