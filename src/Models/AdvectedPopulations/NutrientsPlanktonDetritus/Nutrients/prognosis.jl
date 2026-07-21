
@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:N}, clock, fields, auxiliary_fields) = (
    inorganic_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_nitrogen_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:PO₄}, clock, fields, auxiliary_fields) = (
    inorganic_phosphate_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_phosphate_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields) = (
    inorganic_iron_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_iron_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
)

@inline (bgc::NutrientsNPD)(i, j, k, grid, val_name::Val{:Si}, clock, fields, auxiliary_fields) = (
    inorganic_silicon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_silicon_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
)
