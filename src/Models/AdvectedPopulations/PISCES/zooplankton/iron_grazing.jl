@inline function iron_grazing(g₀,
                               b,
                               pP,
                               pD,
                               pZ,
                               pPOC,
                               J,
                               K,
                               food_threshold_concentration,
                               T,
                               I,
                               food_availability::NamedTuple,
                               iron_availability::NamedTuple)

    base_grazing_rate = g₀ * b ^ T

    total_food = (food_availability.P   * pP +
                  food_availability.D   * pD +
                  food_availability.Z   * pZ +
                  food_availability.POC * pPOC)

    available_total_food = (max(zero(I), food_availability.P   - J) * pP +
                            max(zero(I), food_availability.D   - J) * pD +
                            max(zero(I), food_availability.Z   - J) * pZ +
                            max(zero(I), food_availability.POC - J) * pPOC)

    concentration_limited_grazing = max(zero(I), available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food)

    total_specific_iron_grazing = (
        max(zero(I), food_availability.P   - J) * pP * iron_availability.P +
        max(zero(I), food_availability.D   - J) * pD * iron_availability.D +
        max(zero(I), food_availability.Z   - J) * pZ * iron_availability.Z +
        max(zero(I), food_availability.POC - J) * pPOC * iron_availability.POC
    ) * total_specific_grazing / (available_total_food + eps(zero(I)))

    return total_specific_iron_grazing * I
end

@inline function iron_grazing(zoo::QualityDependantZooplankton, T, I, food_availability::NamedTuple, iron_availability::NamedTuple)
    return iron_grazing(zoo.maximum_grazing_rate,
                        zoo.temperature_sensitivity,
                        zoo.food_preferences.P,
                        zoo.food_preferences.D,
                        zoo.food_preferences.Z,
                        zoo.food_preferences.POC,
                        zoo.specific_food_threshold_concentration,
                        zoo.grazing_half_saturation,
                        zoo.food_threshold_concentration,
                        T,
                        I,
                        food_availability,
                        iron_availability)
end

@inline function iron_grazing(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    food = prey_names(bgc, val_name)
    
    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    food_availability = extract_food_availability(bgc, i, j, k, fields, food)

    iron_ratios = extract_iron_availability(bgc, i, j, k, fields, food)

    return iron_grazing(zoo, T, I, food_availability, iron_ratios) 
end

@inline function iron_flux_feeding(maximum_flux_feeding_rate,
                                   temperature_sensitivity,
                                   T,
                                   I,
                                   sinking_iron_flux)
    g₀ = maximum_flux_feeding_rate
    b = temperature_sensitivity

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_iron_flux

    return total_specific_flux_feeding * I
end

@inline function iron_flux_feeding(zoo::QualityDependantZooplankton, T, I, sinking_iron_flux)
    return iron_flux_feeding(zoo.maximum_flux_feeding_rate,
                             zoo.temperature_sensitivity,
                             T,
                             I,
                             sinking_iron_flux)
end

@inline function iron_flux_feeding(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = edible_iron_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    return iron_flux_feeding(zoo, T, I, sinking_flux)
end
