
@inline function upper_trophic_excretion(zoo, val_name, i, j, k, grid, bgc, clock, fields)
    γ = zoo.dissolved_excretion_fraction

    R = upper_trophic_respiration_product(zoo, val_name, i, j, k, grid, bgc, clock, fields)

    return (1 - γ) * R
end

@inline function upper_trophic_respiration(zoo, val_name, i, j, k, grid, bgc, clock, fields)
    γ = zoo.dissolved_excretion_fraction

    R = upper_trophic_respiration_product(zoo, val_name, i, j, k, grid, bgc, clock, fields)

    return γ * R
end

@inline upper_trophic_respiration_product(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields) = 
    (1 - zoo.minimum_growth_efficiency - zoo.non_assililated_fraction) * upper_trophic_waste(zoo, val_name, i, j, k, grid, bgc, clock, fields)

@inline upper_trophic_fecal_production(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields) =
    zoo.non_assililated_fraction * upper_trophic_waste(zoo, val_name, i, j, k, grid, bgc, clock, fields)

@inline upper_trophic_fecal_iron_production(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields) =
    upper_trophic_fecal_production(zoo, val_name, i, j, k, grid, bgc, clock, fields) * zoo.iron_ratio

@inline function upper_trophic_waste(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    e₀ = zoo.minimum_growth_efficiency
    b  = zoo.temperature_sensetivity
    m₀ = zoo.quadratic_mortality

    temperature_factor = b^T

    I = zooplankton_concentration(val_name, i, j, k, fields)

    return 1 / (1 - e₀) * m₀ * temperature_factor * I^2
end
