using OceanBioME.Models.PISCESModel.Zooplankton: non_assimilated_waste

@inline small_non_assimilated_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    non_assimilated_waste(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields)

@inline large_non_assimilated_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    non_assimilated_waste(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields)

@inline small_non_assimilated_iron_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    non_assimilated_iron_waste(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields)

@inline large_non_assimilated_iron_waste(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    non_assimilated_iron_waste(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields)

@inline small_mortality(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    mortality(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields)

@inline large_mortality(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    linear_mortality(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields)

@inline small_mortality_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    mortality(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields) * zoo.micro.iron_ratio

@inline large_mortality_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields) =
    linear_mortality(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields) * zoo.meso.iron_ratio

@inline total_grazing(zoo::MicroAndMeso, val_prey_name, i, j, k, grid, bgc, clock, fields) =
    (grazing(zoo, val_prey_name, i, j, k, grid, bgc, clock, fields)
     + flux_feeding(zoo, val_prey_name, i, j, k, grid, bgc, clock, fields))
