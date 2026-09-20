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
- `sinking_speeds`: a `NamedTuple` `(sPO = …, bPO = …)` of the small/large particle sinking speeds
  (m/s), used unless `dissolution_lengths` is given
- `dissolution_lengths`: a `NamedTuple` `(sPO = …, bPO = …)` of the small/large particle dissolution
  length scales (m), enabling implicit sinking in place of `sinking_speeds`
- `open_bottom`: whether particulate detritus can sink out of the bottom of the domain
"""
struct CarbonNitrogenDissolvedParticulate{FT, SK} <: AbstractSinkingDetritus{SK}
                       dissolved_remineralisation_rate :: FT
                  small_particle_remineralisation_rate :: FT
                  large_particle_remineralisation_rate :: FT

                         small_fraction_of_solid_waste :: FT

    small_particle_remineralisation_dissolved_fraction :: FT
    large_particle_remineralisation_dissolved_fraction :: FT

                                               sinking :: SK
end

const NPD_CNDP{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate}

required_biogeochemical_tracers(::CarbonNitrogenDissolvedParticulate) =
    (:DON, :DOC, :sPON, :sPOC, :bPON, :bPOC)

required_biogeochemical_tracers(::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}) =
    (:DON, :DOC)

required_biogeochemical_auxiliary_fields(::CarbonNitrogenDissolvedParticulate) = tuple()

function CarbonNitrogenDissolvedParticulate(grid::AbstractGrid{FT};
                                            dissolved_remineralisation_rate = 3.86e-7,
                                            small_particle_remineralisation_rate = 5.88e-7,
                                            large_particle_remineralisation_rate = 5.88e-7,
                                            small_fraction_of_solid_waste = 0.5,
                                            small_particle_remineralisation_dissolved_fraction = 1.0,
                                            large_particle_remineralisation_dissolved_fraction = 1.0,
                                            sinking_speeds = (sPO = 3/day, bPO = 200/day),
                                            dissolution_lengths = nothing,
                                            open_bottom = true) where FT

    if !isnothing(dissolution_lengths)
        dl = (sPON = dissolution_lengths.sPO,
              sPOC = dissolution_lengths.sPO,
              bPON = dissolution_lengths.bPO,
              bPOC = dissolution_lengths.bPO)
        sinking = ImplicitSinking(grid, dl; open_bottom)
    elseif !isnothing(sinking_speeds)
        sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom; three_D = true)
        sinking = ExplicitSinking(sinking_velocities)
    else
        throw(ArgumentError("Must specify either `sinking_speeds` or `dissolution_lengths`"))
    end

    SK = typeof(sinking)

    return CarbonNitrogenDissolvedParticulate{FT, SK}(
        convert(FT, dissolved_remineralisation_rate),
        convert(FT, small_particle_remineralisation_rate),
        convert(FT, large_particle_remineralisation_rate),
        convert(FT, small_fraction_of_solid_waste),
        convert(FT, small_particle_remineralisation_dissolved_fraction),
        convert(FT, large_particle_remineralisation_dissolved_fraction),
        sinking
    )
end

@inline particulate_to_dissolved_nitrogen(i, j, k, d::CarbonNitrogenDissolvedParticulate, fields) = @inbounds (
    fields.sPON[i, j, k] * d.small_particle_remineralisation_rate * d.small_particle_remineralisation_dissolved_fraction +
    fields.bPON[i, j, k] * d.large_particle_remineralisation_rate * d.large_particle_remineralisation_dissolved_fraction
)

@inline particulate_to_dissolved_nitrogen(i, j, k, d::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}, fields) = @inbounds (
    d.sinking.remineralisation.sPON[i, j, k] * d.small_particle_remineralisation_dissolved_fraction +
    d.sinking.remineralisation.bPON[i, j, k] * d.large_particle_remineralisation_dissolved_fraction
)

@inline particulate_to_dissolved_carbon(i, j, k, d::CarbonNitrogenDissolvedParticulate, fields) = @inbounds (
    fields.sPOC[i, j, k] * d.small_particle_remineralisation_rate * d.small_particle_remineralisation_dissolved_fraction +
    fields.bPOC[i, j, k] * d.large_particle_remineralisation_rate * d.large_particle_remineralisation_dissolved_fraction
)

@inline particulate_to_dissolved_carbon(i, j, k, d::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}, fields) = @inbounds (
    d.sinking.remineralisation.sPOC[i, j, k] * d.small_particle_remineralisation_dissolved_fraction +
    d.sinking.remineralisation.bPOC[i, j, k] * d.large_particle_remineralisation_dissolved_fraction
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:DON}, clock, fields, auxiliary_fields) = @inbounds (
    dissolved_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + particulate_to_dissolved_nitrogen(i, j, k, bgc.detritus, fields)
  - grazing(i, j, k, grid, Val(:DON), bgc.plankton, bgc, fields, auxiliary_fields)
  - bgc.detritus.dissolved_remineralisation_rate * fields.DON[i, j, k]
)

@inline (bgc::NPD_CNDP)(i, j, k, grid, ::Val{:DOC}, clock, fields, auxiliary_fields) = @inbounds (
    dissolved_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + particulate_to_dissolved_carbon(i, j, k, bgc.detritus, fields)
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

@inline biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate{<:Any, <:ExplicitSinking}}, ::Union{Val{:sPON}, Val{:sPOC}}) =
    bgc.detritus.sinking.sinking_speeds.sPO

@inline biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate{<:Any, <:ExplicitSinking}}, ::Union{Val{:bPON}, Val{:bPOC}}) =
    bgc.detritus.sinking.sinking_speeds.bPO

@inline inorganic_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DON[i, j, k] * detritus.dissolved_remineralisation_rate
  + fields.sPON[i, j, k] * detritus.small_particle_remineralisation_rate * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + fields.bPON[i, j, k] * detritus.large_particle_remineralisation_rate * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

@inline inorganic_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DON[i, j, k] * detritus.dissolved_remineralisation_rate
  + detritus.sinking.remineralisation.sPON[i, j, k] * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + detritus.sinking.remineralisation.bPON[i, j, k] * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

@inline calcium_carbonate_dissolution(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) * (
        dissolved_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      + fields.sPOC[i, j, k] * detritus.small_particle_remineralisation_rate
      + fields.bPOC[i, j, k] * detritus.large_particle_remineralisation_rate
    )
)

@inline calcium_carbonate_dissolution(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}, bgc, fields, auxiliary_fields) = @inbounds (
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) * (
        dissolved_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      + detritus.sinking.remineralisation.sPOC[i, j, k]
      + detritus.sinking.remineralisation.bPOC[i, j, k]
    )
)

@inline inorganic_carbon_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DOC[i, j, k] * detritus.dissolved_remineralisation_rate
  + fields.sPOC[i, j, k] * detritus.small_particle_remineralisation_rate * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + fields.bPOC[i, j, k] * detritus.large_particle_remineralisation_rate * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
)

@inline inorganic_carbon_waste(i, j, k, grid, detritus::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}, bgc, fields, auxiliary_fields) = @inbounds (
    fields.DOC[i, j, k] * detritus.dissolved_remineralisation_rate
  + detritus.sinking.remineralisation.sPOC[i, j, k] * (1 - detritus.small_particle_remineralisation_dissolved_fraction)
  + detritus.sinking.remineralisation.bPOC[i, j, k] * (1 - detritus.large_particle_remineralisation_dissolved_fraction)
)

@inline implicit_sinking_production(i, j, k, grid, d::CarbonNitrogenDissolvedParticulate, bgc, fields, aux, ::Val{:sPON}) =
    solid_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) * d.small_fraction_of_solid_waste

@inline implicit_sinking_production(i, j, k, grid, d::CarbonNitrogenDissolvedParticulate, bgc, fields, aux, ::Val{:sPOC}) =
    solid_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) * d.small_fraction_of_solid_waste

@inline implicit_sinking_production(i, j, k, grid, d::CarbonNitrogenDissolvedParticulate, bgc, fields, aux, ::Val{:bPON}) =
    solid_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) * (1 - d.small_fraction_of_solid_waste)

@inline implicit_sinking_production(i, j, k, grid, d::CarbonNitrogenDissolvedParticulate, bgc, fields, aux, ::Val{:bPOC}) =
    solid_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) * (1 - d.small_fraction_of_solid_waste)

function Adapt.adapt_structure(to, detritus::CarbonNitrogenDissolvedParticulate{FT}) where FT
    sinking = adapt(to, detritus.sinking)
    SK = typeof(sinking)
    return CarbonNitrogenDissolvedParticulate{FT, SK}(
        detritus.dissolved_remineralisation_rate,
        detritus.small_particle_remineralisation_rate,
        detritus.large_particle_remineralisation_rate,
        detritus.small_fraction_of_solid_waste,
        detritus.small_particle_remineralisation_dissolved_fraction,
        detritus.large_particle_remineralisation_dissolved_fraction,
        sinking
    )
end

Base.summary(::CarbonNitrogenDissolvedParticulate{<:Any, <:ExplicitSinking}) =
    "CarbonNitrogenDissolvedParticulate (DON, DOC, sPON, sPOC, bPON, bPOC)"

Base.summary(::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking}) =
    "CarbonNitrogenDissolvedParticulate (DON, DOC; implicit sPO/bPO sinking)"

function Base.show(io::IO, dp::CarbonNitrogenDissolvedParticulate{<:Any, <:ExplicitSinking})
    msg  = "CarbonNitrogenDissolvedParticulate\n"
    msg *= "└── Particle sinking speeds\n"
    msg *= "  ├── sPOX : " * summary(dp.sinking.sinking_speeds.sPO.w) * "\n"
    msg *= "  └── bPOX : " * summary(dp.sinking.sinking_speeds.bPO.w)
    print(io, msg)
    return nothing
end

function Base.show(io::IO, dp::CarbonNitrogenDissolvedParticulate{<:Any, <:ImplicitSinking})
    msg  = "CarbonNitrogenDissolvedParticulate\n"
    msg *= "└── Sinking: $(summary(dp.sinking))"
    print(io, msg)
    return nothing
end
