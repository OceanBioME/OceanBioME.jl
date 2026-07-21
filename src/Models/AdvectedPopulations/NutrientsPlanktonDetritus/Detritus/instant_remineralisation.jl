"""
    InstantRemineralisationDetritus()

The simplest detritus component for the `detritus` slot of a [`NutrientsPlanktonDetritus`](@ref)
model. It adds no tracers of its own: plankton dissolved and solid waste is instantaneously
remineralised straight back to the inorganic nutrient pool (and, when present, inorganic carbon and
calcite). This is the default detritus and is useful for closed-budget box models or when organic
matter export is not of interest.
"""
struct InstantRemineralisationDetritus end

required_biogeochemical_tracers(::InstantRemineralisationDetritus) = tuple()
required_biogeochemical_auxiliary_fields(::InstantRemineralisationDetritus) = tuple()

@inline inorganic_waste(i, j, k, grid, ::InstantRemineralisationDetritus, bgc, args...) =
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, args...) + solid_waste(i, j, k, grid, bgc.plankton, bgc, args...)

@inline calcite_dissolution(i, j, k, grid, ::InstantRemineralisationDetritus, bgc, fields, auxiliary_fields) = (
    inorganic_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)
