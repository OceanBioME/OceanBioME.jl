"""
This module adds active particles which sink, dissolve, and release material into tracer fields. For example, mineral particles added to the ocean for the purpose of ocean alkalinity enhancement (OAE) or dissolving secondarary mineral precipitate.

The particles are defined by a radius, r (m), which sets how fast they dissolve and sink.

Prognostic properties
=====================
* Radius: r (m)

Tracer dependencies
===================
None

Coupled tracers
===============
Those listed in `released_tracers`, by default alkalinity: Alk (mmol/m³)
"""

module DissolvingParticlesModel

export DissolvingMineral, DissolvingParticles, StokesSinking

using Oceananigans.Grids: AbstractGrid

using OceanBioME.Particles: BiogeochemicalParticles

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers, update_particle_state!
import Oceananigans.Models.LagrangianParticleTracking: particle_w_velocity

import Base: summary, show

"""
    StokesSinking(FT = Float64;
                  particle_density,
                  water_density = 1025,
                  dynamic_viscosity = 1.08e-3,
                  gravitational_acceleration = 9.81)

Function to calculate Stokes' law for the speed (m/s) at which a sphere of radius `r` (m) sinks:

```math
w = \\frac{2(\\rho_p - \\rho)gr^2}{9\\mu},
```

where ``\\rho_p`` is the `particle_density` (kg/m³), ``\\rho`` the `water_density` (kg/m³),
``\\mu`` the `dynamic_viscosity` (kg/m/s), and ``g`` the `gravitational_acceleration` (m/s²).

This is only valid when the particle Reynolds number, ``2rw\\rho/\\mu``, is small,
i.e. for particles smaller than about 100 μm.
"""
struct StokesSinking{FT}
              particle_density :: FT
                 water_density :: FT
             dynamic_viscosity :: FT
    gravitational_acceleration :: FT
end

StokesSinking(FT = Float64;
              particle_density,
              water_density = 1025,
              dynamic_viscosity = 1.08e-3,
              gravitational_acceleration = 9.81) =
    StokesSinking{FT}(particle_density, water_density, dynamic_viscosity, gravitational_acceleration)

@inline function (stokes::StokesSinking)(r)
    ρₚ = stokes.particle_density
    ρ  = stokes.water_density
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
                      sinking_speed = StokesSinking(FT; particle_density = density))

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

Keyword Arguments
=================

- `density`: density of the mineral (kg/m³)
- `molar_mass`: molar mass of the mineral (kg/mol)
- `released_tracers`: `NamedTuple` of the moles of each tracer released per mole of mineral dissolved,
  e.g. `(; Alk = 2, DIC = 1)` for calcium carbonate
- `dissolution_rate`: the speed (m/s) that the radius shrinks, ``k(r)``, either a number or a function of radius
- `sinking_speed`: the speed (m/s, positive downwards) that the particles sink, either a number or a function of radius

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
└── sinking speed: StokesSinking{Float64}(2710.0, 1025.0, 0.00108, 9.81) m/s
```
"""
struct DissolvingMineral{FT, R, D, W}
             density :: FT
          molar_mass :: FT
    released_tracers :: R
    dissolution_rate :: D
       sinking_speed :: W
end

function DissolvingMineral(FT = Float64;
                           density = 2211,
                           molar_mass = 74.093e-3,
                           released_tracers = (; Alk = 2),
                           dissolution_rate = 1e-8,
                           sinking_speed = StokesSinking(FT; particle_density = density))

    released_tracers = map(ν -> convert(FT, ν), released_tracers)
    dissolution_rate = maybe_convert(FT, dissolution_rate)
    sinking_speed = maybe_convert(FT, sinking_speed)

    return DissolvingMineral(convert(FT, density), convert(FT, molar_mass), released_tracers, dissolution_rate, sinking_speed)
end

maybe_convert(FT, value::Number) = convert(FT, value)
maybe_convert(FT, func) = func

# `dissolution_rate` and `sinking_speed` can either be numbers or functions of radius
@inline evaluate(value::Number, r) = value
@inline evaluate(func, r) = func(r)

"""
    DissolvingParticles(n; grid, mineral_parameters = NamedTuple(), kwargs...)

Sets up `n` dissolving `BiogeochemicalParticles` with default [`DissolvingMineral`](@ref) parameters except those specified in `mineral_parameters`. `kwargs` are passed onto `BiogeochemicalParticles`.
"""
DissolvingParticles(n; grid::AbstractGrid{FT}, mineral_parameters = NamedTuple(), kwargs...) where FT =
    BiogeochemicalParticles(n; grid, biogeochemistry = DissolvingMineral(FT; mineral_parameters...), kwargs...)

@inline required_particle_fields(::DissolvingMineral) = (:r, )
@inline required_tracers(::DissolvingMineral) = ()
@inline coupled_tracers(mineral::DissolvingMineral) = keys(mineral.released_tracers)

# Volume of mineral dissolving per second (m³/s), which stops once the particle has fully dissolved
@inline dissolution(mineral, r) = ifelse(r > 0, 4π * r^2 * evaluate(mineral.dissolution_rate, r), zero(r))

# Change in radius (m/s)
@inline (mineral::DissolvingMineral)(::Val{:r}, t, r) =
    - ifelse(r > 0, evaluate(mineral.dissolution_rate, r), zero(r))

# Release of each tracer (mmol/s)
@inline function (mineral::DissolvingMineral)(::Val{tracer}, t, r) where tracer
    molar_volume = mineral.molar_mass / mineral.density # m³/mol
    ν = mineral.released_tracers[tracer] # mol tracer / mol mineral

    return 1000 * ν * dissolution(mineral, r) / molar_volume
end

@inline sinking_speed(mineral, r) = ifelse(r > 0, evaluate(mineral.sinking_speed, r), zero(r))

# Oceananigans' particle advection asks for each particle's vertical velocity,
# so the particles sink relative to the water when they're moved by `LagrangianAdvection`
@inline particle_w_velocity(particles::BiogeochemicalParticles{N, <:DissolvingMineral}, p, w) where N =
    w - sinking_speed(particles.biogeochemistry, @inbounds(particles.fields.r[p]))

# The final timestep of dissolution can overshoot to a negative radius, so we clip it
# (this doesn't change what was released because nothing is released when `r ≤ 0`)
function update_particle_state!(particles::BiogeochemicalParticles{N, <:DissolvingMineral}, model, Δt) where N
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
              "└── sinking speed: ", mineral.sinking_speed, " m/s")

end # module
