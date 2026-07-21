using Oceananigans.Units

using ..NutrientsPlanktonDetritusModels: 
    dissolved_nitrogen_waste, 
    dissolved_carbon_waste, 
    solid_nitrogen_waste, 
    solid_carbon_waste

import ..NutrientsPlanktonDetritusModels.InorganicCarbonModels: 
    inorganic_carbon_waste

"""
    CarbonNitrogenDissolvedParticulate(grid; kwargs...)

A variable-Redfield detritus component for the `detritus` slot of a
[`NutrientsPlanktonDetritus`](@ref) model that tracks carbon and nitrogen separately in dissolved,
small-particulate, and large-particulate classes. It adds the tracers `DON`, `DOC`, `sPON`, `sPOC`,
`bPON`, and `bPOC`; the two particulate classes sink, and each pool remineralises to inorganic
nutrients and carbon (partly via the dissolved pool).

Keyword Arguments
=================

- `grid`: (required) the geometry, needed to configure the sinking-speed fields
- `dissolved_remineralisation_rate`, `small_particle_remineralisation_rate`,
  `large_particle_remineralisation_rate`: per-class remineralisation rates (1/s)
- `small_fraction_of_solid_waste`: the fraction of solid plankton waste routed to the small
  particulate class (the rest goes to the large class)
- `small_particle_remineralisation_dissolved_fraction`,
  `large_particle_remineralisation_dissolved_fraction`: the fraction of each particulate class's
  remineralisation that passes through the dissolved pool
- `sinking_speeds`: a `NamedTuple` `(sPO = …, bPO = …)` of the small/large particle sinking speeds (m/s)
- `open_bottom`: whether particulate detritus can sink out of the bottom of the domain
"""
struct CarbonNitrogenDissolvedParticulate{FT, SV}
                       dissolved_remineralisation_rate :: FT
                  small_particle_remineralisation_rate :: FT
                  large_particle_remineralisation_rate :: FT

                         small_fraction_of_solid_waste :: FT

    small_particle_remineralisation_dissolved_fraction :: FT
    large_particle_remineralisation_dissolved_fraction :: FT

                                    sinking_velocities :: SV
end

const NPD_CNDP{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate}

required_biogeochemical_tracers(::CarbonNitrogenDissolvedParticulate) =
    (:DON, :DOC, :sPON, :sPOC, :bPON, :bPOC)

required_biogeochemical_auxiliary_fields(::CarbonNitrogenDissolvedParticulate) = tuple()

function CarbonNitrogenDissolvedParticulate(grid::AbstractGrid{FT};
                                            dissolved_remineralisation_rate = 3.86e-7,
                                            small_particle_remineralisation_rate = 5.88e-7,
                                            large_particle_remineralisation_rate = 5.88e-7,
                                            small_fraction_of_solid_waste = 0.5,
                                            small_particle_remineralisation_dissolved_fraction = 1.0,
                                            large_particle_remineralisation_dissolved_fraction = 1.0,
                                            sinking_speeds = (sPO = 3/day, bPO = 200/day),
                                            open_bottom = true) where FT

    sinking_velocities = setup_velocity_fields((; sPO = convert(FT, sinking_speeds.sPO),
                                                  bPO = convert(FT, sinking_speeds.bPO)),
                                               grid, open_bottom; three_D = true)
    SV = typeof(sinking_velocities)

    return CarbonNitrogenDissolvedParticulate{FT, SV}(
        convert(FT, dissolved_remineralisation_rate),
        convert(FT, small_particle_remineralisation_rate),
        convert(FT, large_particle_remineralisation_rate),
        convert(FT, small_fraction_of_solid_waste),
        convert(FT, small_particle_remineralisation_dissolved_fraction),
        convert(FT, large_particle_remineralisation_dissolved_fraction),
        sinking_velocities
    )
end

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:DON}, clock, fields, auxiliary_fields) = @inbounds (
    dissolved_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + fields.sPON[i, j, k] * bgc.detritus.small_particle_remineralisation_rate * 
      bgc.detritus.small_particle_remineralisation_dissolved_fraction
  + fields.bPON[i, j, k] * bgc.detritus.large_particle_remineralisation_rate * 
      bgc.detritus.large_particle_remineralisation_dissolved_fraction
  - grazing(i, j, k, grid, Val(:DON), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.dissolved_remineralisation_rate * fields.DON[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:DOC}, clock, fields, auxiliary_fields) = @inbounds (
    dissolved_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + fields.sPOC[i, j, k] * bgc.detritus.small_particle_remineralisation_rate * 
        bgc.detritus.small_particle_remineralisation_dissolved_fraction
  + fields.bPOC[i, j, k] * bgc.detritus.large_particle_remineralisation_rate * 
        bgc.detritus.large_particle_remineralisation_dissolved_fraction
  - grazing(i, j, k, grid, Val(:DOC), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.dissolved_remineralisation_rate * fields.DOC[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:sPON}, clock, fields, auxiliary_fields) = @inbounds (
    solid_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * bgc.detritus.small_fraction_of_solid_waste
  - grazing(i, j, k, grid, Val(:sPON), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.small_particle_remineralisation_rate * fields.sPON[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:sPOC}, clock, fields, auxiliary_fields) = @inbounds (
    solid_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * bgc.detritus.small_fraction_of_solid_waste
  - grazing(i, j, k, grid, Val(:sPOC), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.small_particle_remineralisation_rate * fields.sPOC[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:bPON}, clock, fields, auxiliary_fields) = @inbounds (
    solid_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * (1 - bgc.detritus.small_fraction_of_solid_waste)
  - grazing(i, j, k, grid, Val(:bPON), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.large_particle_remineralisation_rate * fields.bPON[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:bPOC}, clock, fields, auxiliary_fields) = @inbounds (
    solid_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * (1 - bgc.detritus.small_fraction_of_solid_waste)
  - grazing(i, j, k, grid, Val(:bPOC), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.large_particle_remineralisation_rate * fields.bPOC[i, j, k]
)

@inline biogeochemical_drift_velocity(bgc::NPD_CNDP, ::Union{Val{:sPON}, Val{:sPOC}}) = bgc.detritus.sinking_velocities.sPO
@inline biogeochemical_drift_velocity(bgc::NPD_CNDP, ::Union{Val{:bPON}, Val{:bPOC}}) = bgc.detritus.sinking_velocities.bPO

@inline inorganic_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DON[i, j, k] * detritus.dissolved_remineralisation_rate
  + fields.sPON[i, j, k] * detritus.small_particle_remineralisation_rate * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + fields.bPON[i, j, k] * detritus.large_particle_remineralisation_rate * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

@inline calcite_dissolution(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) * (
        dissolved_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      + fields.sPOC[i, j, k] * detritus.small_particle_remineralisation_rate
      + fields.bPOC[i, j, k] * detritus.large_particle_remineralisation_rate
    )
)

@inline inorganic_carbon_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DOC[i, j, k] * detritus.dissolved_remineralisation_rate
  + fields.sPOC[i, j, k] * detritus.small_particle_remineralisation_rate * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + fields.bPOC[i, j, k] * detritus.large_particle_remineralisation_rate * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
)

function Adapt.adapt_structure(to, detritus::CarbonNitrogenDissolvedParticulate{FT}) where FT
    sinking_velocities = adapt(to, detritus.sinking_velocities)
    SV = typeof(sinking_velocities)
    return CarbonNitrogenDissolvedParticulate{FT, SV}(
        detritus.dissolved_remineralisation_rate,
        detritus.small_particle_remineralisation_rate,
        detritus.large_particle_remineralisation_rate,
        detritus.small_fraction_of_solid_waste,
        detritus.small_particle_remineralisation_dissolved_fraction,
        detritus.large_particle_remineralisation_dissolved_fraction,
        sinking_velocities
    )
end

Base.summary(::CarbonNitrogenDissolvedParticulate) =
    "CarbonNitrogenDissolvedParticulate (DON, DOC, sPON, sPOC, bPON, bPOC)"

function Base.show(io::IO, dp::CarbonNitrogenDissolvedParticulate)
    msg  = "CarbonNitrogenDissolvedParticulate\n"
    msg *= "└── Particle sinking speeds\n"
    msg *= "  ├── sPOX : " * summary(dp.sinking_velocities.sPO.w) * "\n"
    msg *= "  └── bPOX : " * summary(dp.sinking_velocities.bPO.w)
    print(io, msg)
    return nothing
end
