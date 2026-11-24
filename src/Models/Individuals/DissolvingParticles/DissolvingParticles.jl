"""
Active particles which sink, dissolve, and release material into tracer fields, for example mineral
particles added to the ocean for ocean alkalinity enhancement, or dissolving secondary mineral precipitate.

Prognostic properties
=====================
* Radius: r (m)

Diagnostic properties
=====================
* Vertical velocity relative to the water: sinking_speed (m/s, positive upwards)

Tracer dependencies
===================
Those listed in `required_fields`, by default none

Coupled tracers
===============
Those listed in `released_tracers`, by default alkalinity: Alk (mmol/m³)
"""

module DissolvingParticlesModel

export DissolvingMineral, DissolvingParticles, StokesSinking

using Oceananigans: Biogeochemistry, fields
using Oceananigans.Architectures: architecture
using Oceananigans.Models: LagrangianParticleTracking
using Oceananigans.Utils: KernelParameters, launch!

using OceanBioME: Particles
using OceanBioME.Particles: BiogeochemicalParticles, extract_tracer_values, required_tracers

using KernelAbstractions: @kernel, @index

"""
    StokesSinking(FT = Float64;
                  particle_density,
                  water_density = 1025,
                  dynamic_viscosity = 1.08e-3,
                  gravitational_acceleration = 9.81)

Stokes' law for the vertical velocity (m/s, positive upwards) of a sphere of radius `r` (m) relative to the water,

    w = 2 (ρ - ρₚ) g r² / 9μ,

valid at small particle Reynolds number, `2 r |w| ρ / μ`. It is called as
`stokes(mineral, particle_fields, tracer_fields)` and reads the radius from `particle_fields.r`.

Keyword Arguments
=================

- `particle_density`: density of the particle, `ρₚ` (kg/m³).
- `water_density`: density of the water, `ρ` (kg/m³), a number or a function called as
  `water_density(mineral, particle_fields, tracer_fields)`. Default: 1025.
- `dynamic_viscosity`: dynamic viscosity of the water, `μ` (kg/m/s), a number or a function called as
  `dynamic_viscosity(mineral, particle_fields, tracer_fields)`. Default: 1.08e-3.
- `gravitational_acceleration`: `g` (m/s²). Default: 9.81.

Example
=======

```jldoctest
julia> using OceanBioME

julia> stokes = StokesSinking(particle_density = 2710)
StokesSinking{Float64, Float64, Float64}(2710.0, 1025.0, 0.00108, 9.81)
```
"""
struct StokesSinking{WD, DV, FT}
              particle_density :: FT
                 water_density :: WD
             dynamic_viscosity :: DV
    gravitational_acceleration :: FT
end

function StokesSinking(FT = Float64;
                       particle_density,
                       water_density = 1025,
                       dynamic_viscosity = 1.08e-3,
                       gravitational_acceleration = 9.81)

    water_density = maybe_convert(FT, water_density)
    dynamic_viscosity = maybe_convert(FT, dynamic_viscosity)

    return StokesSinking(convert(FT, particle_density), water_density, dynamic_viscosity,
                         convert(FT, gravitational_acceleration))
end

@inline function (stokes::StokesSinking)(mineral, particle_fields, tracer_fields)
    r  = particle_fields.r
    ρₚ = stokes.particle_density
    ρ  = evaluate(stokes.water_density, mineral, particle_fields, tracer_fields)
    μ  = evaluate(stokes.dynamic_viscosity, mineral, particle_fields, tracer_fields)
    g  = stokes.gravitational_acceleration

    return 2 * (ρ - ρₚ) * g * r^2 / (9 * μ)
end

"""
    DissolvingMineral(FT = Float64;
                      density = 2211,
                      molar_mass = 74.093e-3,
                      released_tracers = (; DIC = 1, Alk = 2),
                      dissolution_rate = 1e-8,
                      sinking_speed = StokesSinking(FT; particle_density = density),
                      required_fields = ())

Spherical mineral particles of radius `r` (m) which dissolve, release tracers, and sink.

The radius shrinks at the dissolution rate `k` (m/s),

    dr/dt = -k,

so each particle releases `4π r² k θ ρₚ / M` of each tracer (mol/s, converted to mmol/s),
where `θ` is the moles of tracer released per mole of mineral, `ρₚ` the `density`,
and `M` the `molar_mass`.

`dissolution_rate` and `sinking_speed` are numbers or functions called as
`f(mineral, particle_fields, tracer_fields)`, where `particle_fields` is a `NamedTuple` of the particle's
`r` and `sinking_speed`, and `tracer_fields` is a `NamedTuple` of the `required_fields` at the particle position.

Keyword Arguments
=================

- `density`: density of the mineral (kg/m³). Default: 2211.
- `molar_mass`: molar mass of the mineral (kg/mol). Default: 74.093e-3.
- `released_tracers`: `NamedTuple` of the moles of each tracer released per mole of mineral dissolved.
  Default: `(; DIC = 1, Alk = 2)`.
- `dissolution_rate`: speed that the radius shrinks (m/s). Default: 1e-8.
- `sinking_speed`: vertical velocity of the particles relative to the water (m/s, positive upwards).
  Default: `StokesSinking(FT; particle_density = density)`.
- `required_fields`: names of the model fields (tracers, velocities, or auxiliary fields) passed to
  `dissolution_rate` and `sinking_speed`. Default: `()`.

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
├── sinking speed: StokesSinking{Float64, Float64, Float64}(2710.0, 1025.0, 0.00108, 9.81) m/s
└── required fields: ()
```
"""
struct DissolvingMineral{FT, R, D, W, F}
             density :: FT
          molar_mass :: FT
    released_tracers :: R
    dissolution_rate :: D
       sinking_speed :: W

    DissolvingMineral(density::FT, molar_mass::FT, released_tracers::R, dissolution_rate::D, sinking_speed::W,
                      required_fields::Tuple) where {FT, R, D, W} =
        new{FT, R, D, W, required_fields}(density, molar_mass, released_tracers, dissolution_rate, sinking_speed)
end

function DissolvingMineral(FT = Float64;
                           density = 2211,
                           molar_mass = 74.093e-3,
                           released_tracers = (; DIC = 1, Alk = 2),
                           dissolution_rate = 1e-8,
                           sinking_speed = StokesSinking(FT; particle_density = density),
                           required_fields = ())

    released_tracers = map(θ -> convert(FT, θ), released_tracers)
    dissolution_rate = maybe_convert(FT, dissolution_rate)
    sinking_speed = maybe_convert(FT, sinking_speed)

    return DissolvingMineral(convert(FT, density), convert(FT, molar_mass), released_tracers,
                             dissolution_rate, sinking_speed, required_fields)
end

maybe_convert(FT, value::Number) = convert(FT, value)
maybe_convert(FT, func) = func

@inline evaluate(value::Number, args...) = value
@inline evaluate(func, args...) = func(args...)

"""
    DissolvingParticles(n;
                        grid,
                        mineral_parameters = NamedTuple(),
                        biogeochemistry = DissolvingMineral(eltype(grid); mineral_parameters...),
                        kwargs...)

Return `n` `BiogeochemicalParticles` on `grid` with `DissolvingMineral` `biogeochemistry`.
`mineral_parameters` are keyword arguments for `DissolvingMineral`, and `kwargs` are passed to
`BiogeochemicalParticles`.
"""
DissolvingParticles(n;
                    grid,
                    mineral_parameters = NamedTuple(),
                    biogeochemistry = DissolvingMineral(eltype(grid); mineral_parameters...),
                    kwargs...) =
    BiogeochemicalParticles(n; grid, biogeochemistry, kwargs...)

@inline Particles.required_particle_fields(::DissolvingMineral) = (:r, :sinking_speed)
@inline Particles.required_tracers(::DissolvingMineral{<:Any, <:Any, <:Any, <:Any, F}) where F = F
@inline Particles.coupled_tracers(mineral::DissolvingMineral) = keys(mineral.released_tracers)

@inline function dissolution_rate(mineral, r, w, field_values...)
    particle_fields = (; r, sinking_speed = w)
    tracer_fields = NamedTuple{required_tracers(mineral)}(field_values)

    return evaluate(mineral.dissolution_rate, mineral, particle_fields, tracer_fields) * (r > 0)
end

@inline (mineral::DissolvingMineral)(::Val{:r}, t, r, w, field_values...) =
    - dissolution_rate(mineral, r, w, field_values...)

@inline (mineral::DissolvingMineral)(::Val{:sinking_speed}, t, r, w, field_values...) = zero(r)

@inline dissolution(mineral, r::FT, w, field_values...) where FT =
    convert(FT, 4π) * r^2 * dissolution_rate(mineral, r, w, field_values...) # m³ / s

@inline function (mineral::DissolvingMineral{FT})(::Val{tracer}, t, r, w, field_values...) where {FT, tracer}
    mineral_concentration = convert(FT, 1000) * mineral.density / mineral.molar_mass # mmol / m³
    θ = mineral.released_tracers[tracer] # mol tracer / mol mineral
    dissolution_rate = dissolution(mineral, r, w, field_values...) # m³ / s

    return θ * dissolution_rate * mineral_concentration # mmol / s
end

@inline sinking_speed(mineral, particle_fields, tracer_fields) =
    evaluate(mineral.sinking_speed, mineral, particle_fields, tracer_fields) * (particle_fields.r > 0)

const DissolvingBiogeochemicalParticles = BiogeochemicalParticles{<:Any, <:DissolvingMineral}

@inline LagrangianParticleTracking.particle_w_velocity(particles::DissolvingBiogeochemicalParticles, n, w) =
    @inbounds w + particles.fields.sinking_speed[n]

function Particles.update_particle_state!(particles::DissolvingBiogeochemicalParticles, model, Δt)
    r = particles.fields.r
    r .= max.(r, 0)

    return nothing
end

Biogeochemistry.update_biogeochemical_state!(model, particles::DissolvingBiogeochemicalParticles) =
    update_sinking_speed!(particles, model)

function update_sinking_speed!(particles, model)
    grid = model.grid

    launch!(architecture(grid), grid, KernelParameters(1:length(particles)),
            _update_sinking_speed!, particles, grid, fields(model))

    return nothing
end

@kernel function _update_sinking_speed!(particles, grid, fields)
    n = @index(Global)

    mineral = particles.biogeochemistry

    particle_fields = map(field -> @inbounds(field[n]), particles.fields)

    tracer_values = extract_tracer_values(nothing, particles.field_interpolation, particles, grid, fields, n)
    tracer_fields = NamedTuple{required_tracers(mineral)}(tracer_values)

    @inbounds particles.fields.sinking_speed[n] = sinking_speed(mineral, particle_fields, tracer_fields)
end

Base.summary(::DissolvingMineral{FT}) where FT = "DissolvingMineral{$FT} biogeochemistry"

Base.show(io::IO, mineral::DissolvingMineral) =
    print(io, summary(mineral), " tracking particle `r`adius", "\n",
              "├── density: ", mineral.density, " kg/m³", "\n",
              "├── molar mass: ", mineral.molar_mass, " kg/mol", "\n",
              "├── released tracers: ", mineral.released_tracers, " mol/mol", "\n",
              "├── dissolution rate: ", mineral.dissolution_rate, " m/s", "\n",
              "├── sinking speed: ", mineral.sinking_speed, " m/s", "\n",
              "└── required fields: ", required_tracers(mineral))

end # module
