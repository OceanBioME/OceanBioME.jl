struct InstantRemineralisation end

required_biogeochemical_tracers(::InstantRemineralisation) = tuple()
required_biogeochemical_auxiliary_fields(::InstantRemineralisation) = tuple()

@inline inorganic_waste(i, j, k, grid, ::InstantRemineralisation, bgc, args...) =
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, args...) + solid_waste(i, j, k, grid, bgc.plankton, bgc, args...)

@inline calcite_dissolution(i, j, k, grid, ::InstantRemineralisation, bgc, fields, auxiliary_fields) = (
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)
