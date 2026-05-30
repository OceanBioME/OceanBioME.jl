
"""
    nutrient_uptake

Generic nutrient uptake expecting arguments
`plankton, bgc, i, j, k, val_name, fields, auxiliary_fields`
"""
function nutrient_uptake end

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:N}, clock, fields, auxiliary_fields) = (
    inorganic_nitrogen_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_nitrogen_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:PO₄}, clock, fields, auxiliary_fields) = (
    inorganic_phosphate_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_phosphate_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields) = (
    inorganic_iron_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_iron_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:Si}, clock, fields, auxiliary_fields) = (
    inorganic_silicon_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_silicon_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
)
