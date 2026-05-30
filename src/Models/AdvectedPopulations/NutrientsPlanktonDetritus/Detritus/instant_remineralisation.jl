struct InstantRemineralisation end

required_biogeochemical_tracers(::InstantRemineralisation) = tuple()
required_biogeochemical_auxiliary_fields(::InstantRemineralisation) = tuple()

@inline inorganic_waste(::InstantRemineralisation, bgc, args...) =
    dissolved_waste(bgc.plankton, bgc, args...) + solid_waste(bgc.plankton, bgc, args...)

@inline calcite_dissolution(::InstantRemineralisation, bgc, i, j, k, fields, auxiliary_fields) = (
    dissolved_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  * carbon_ratio(bgc.plankton, bgc, i, j, k, fields)
  * calcite_rain_ratio(bgc.plankton, bgc, i, j, k, fields)
)
