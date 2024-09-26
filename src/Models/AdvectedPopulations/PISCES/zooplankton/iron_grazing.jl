
@inline function iron_grazing(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    # food quantity
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensetivity
    p    = zoo.food_preferences
    food = prey_names(bgc, val_name)
    J    = zoo.specific_food_thresehold_concentration
    K    = zoo.grazing_half_saturation
    food_threshold_concentration = zoo.food_threshold_concentration

    N = length(food)

    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    base_grazing_rate = g₀ * b ^ T

    food_availability = extract_food_availability(bgc, i, j, k, fields, food)

    total_food = sum(ntuple(n->food_availability[n] * p[n], Val(N)))

    available_total_food = sum(ntuple(n->max(zero(grid), (food_availability[n] - J)) * p[n], Val(N)))

    concentration_limited_grazing = max(0, available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food)

    iron_ratios = extract_iron_availability(bgc, i, j, k, fields, food)

    total_specific_iron_grazing = sum(ntuple(n->max(zero(grid), (food_availability[n] - J)) * p[n] * iron_ratios[n], Val(N))) * total_specific_grazing / (available_total_food + eps(0.0))

    return total_specific_iron_grazing * I
end

@inline function iron_flux_feeding(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀ = zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensetivity

    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = edible_iron_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_flux 

    return total_specific_flux_feeding * I
end

