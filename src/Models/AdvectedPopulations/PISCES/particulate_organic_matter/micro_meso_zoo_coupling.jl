using OceanBioME.Models.PISCESModel.Zooplankton: non_assimilated_waste

@inline small_non_assimilated_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    non_assimilated_waste(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline large_non_assimilated_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    non_assimilated_waste(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline small_non_assimilated_iron_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    non_assimilated_iron_waste(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline large_non_assimilated_iron_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    non_assimilated_iron_waste(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline small_mortality(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    mortality(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline large_mortality(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    linear_mortality(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline small_mortality_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    mortality(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields) * zoo.micro.iron_ratio

@inline large_mortality_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    linear_mortality(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields) * zoo.meso.iron_ratio

@inline total_grazing(zoo::MicroAndMeso, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (grazing(zoo, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + flux_feeding(zoo, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline total_grazing(zoo::MicroAndMeso, val_prey_name::Val{:GOC}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    flux_feeding(zoo, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
