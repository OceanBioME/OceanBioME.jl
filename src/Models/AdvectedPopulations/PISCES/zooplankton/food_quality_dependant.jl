"""
    QualityDependantZooplankton

The PISCES zooplankton growth model where each class has preferences
for grazing on nanophytoplankton (P), diatoms (D), microzooplankton (Z),
and particulate organic matter (POC), and can flux feed on sinking 
particulates (POC and GOC).

This model assumes a fixed ratio for all other elements (i.e. N, P, Fe).
"""
struct QualityDependantZooplankton{FT, FP}
                   temperature_sensitivity :: FT #
                      maximum_grazing_rate :: FT # 1 / s

                          food_preferences :: FP 

              food_threshold_concentration :: FT # mmol C / m³
    specific_food_threshold_concentration :: FT # mmol C / m³

                   grazing_half_saturation :: FT # mmol C / m³

                 maximum_flux_feeding_rate :: FT # m / (mmol C / m³)

                                iron_ratio :: FT # μmol Fe / mmol C

                 minimum_growth_efficiency :: FT #
                  non_assimilated_fraction :: FT #

                 mortality_half_saturation :: FT # mmol C / m³
                       quadratic_mortality :: FT # 1 / (mmol C / m³) / s
                          linear_mortality :: FT # 1 / s

    # this should be called inorganic excretion factor
              dissolved_excretion_fraction :: FT #
              undissolved_calcite_fraction :: FT #

    QualityDependantZooplankton{FT, FP}(temperature_sensitivity, 
                                        maximum_grazing_rate,
                                        food_preferences,
                                        food_threshold_concentration, 
                                        specific_food_threshold_concentration,
                                        grazing_half_saturation, 
                                        maximum_flux_feeding_rate,
                                        iron_ratio,
                                        minimum_growth_efficiency, 
                                        non_assimilated_fraction,
                                        mortality_half_saturation, 
                                        quadratic_mortality, 
                                        linear_mortality,
                                        dissolved_excretion_fraction, 
                                        undissolved_calcite_fraction) where {FT, FP} = 
        new{FT, FP}(temperature_sensitivity, 
                    maximum_grazing_rate,
                    food_preferences,
                    food_threshold_concentration, 
                    specific_food_threshold_concentration,
                    grazing_half_saturation, 
                    maximum_flux_feeding_rate,
                    iron_ratio,
                    minimum_growth_efficiency, 
                    non_assimilated_fraction,
                    mortality_half_saturation, 
                    quadratic_mortality, 
                    linear_mortality,
                    dissolved_excretion_fraction, 
                    undissolved_calcite_fraction)

    function QualityDependantZooplankton(FT = Float64;
                                         temperature_sensitivity = 1.079, #
                                         maximum_grazing_rate, # 1 / s

                                         food_preferences, 

                                         food_threshold_concentration = 0.3, # mmol C / m³
                                         specific_food_threshold_concentration = 0.001, # mmol C / m³

                                         grazing_half_saturation = 20.0, # mmol C / m³

                                         maximum_flux_feeding_rate, # m / (mmol C / m³)

                                         iron_ratio, # μmol Fe / mmol C

                                         minimum_growth_efficiency, #
                                         non_assimilated_fraction = 0.3, #

                                         mortality_half_saturation = 0.2, # mmol C / m³
                                         quadratic_mortality, # 1 / (mmol C / m³) / s
                                         linear_mortality, # 1 / s

                                         # this should be called inorganic excretion factor
                                         dissolved_excretion_fraction = 0.6, #
                                         undissolved_calcite_fraction) #
        
        food_preferences = NamedTuple{keys(food_preferences)}(map(fp -> convert(FT, fp), food_preferences))

        FP = typeof(food_preferences)

        return new{FT, FP}(temperature_sensitivity, maximum_grazing_rate,
                           food_preferences,
                           food_threshold_concentration, specific_food_threshold_concentration,
                           grazing_half_saturation, maximum_flux_feeding_rate,
                           iron_ratio, 
                           minimum_growth_efficiency, non_assimilated_fraction,
                           mortality_half_saturation, quadratic_mortality, linear_mortality,
                           dissolved_excretion_fraction, undissolved_calcite_fraction)
    end
end

required_biogeochemical_tracers(::QualityDependantZooplankton, name_base) = tuple(name_base)

@inline function growth_death(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    gI, e = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 
    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 
    mI = mortality(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return e * (gI + gfI) - mI
end

# fallback
@inline extract_food_availability(bgc, i, j, k, fields, names::NTuple{N}) where N =
    ntuple(n -> concentration(Val(names[n]), i, j, k, fields), Val(N))

@inline extract_iron_availability(bgc, i, j, k, fields, names::NTuple{N}) where N =
    ntuple(n -> iron_ratio(Val(names[n]), i, j, k, bgc, fields), Val(N))

@inline function grazing(g₀,
                                 b,
                                 pP,
                                 pD,
                                 pZ,
                                 pPOC,
                                 J,
                                 K,
                                 food_threshold_concentration,
                                 θFe,
                                 e₀,
                                 σ,
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

    total_iron = (iron_availability.P   * pP +
                  iron_availability.D   * pD +
                  iron_availability.Z   * pZ +
                  iron_availability.POC * pPOC)

    iron_grazing_ratio = total_iron / (θFe * total_specific_grazing + eps(zero(I)))

    food_quality = min(one(I), iron_grazing_ratio)

    growth_efficiency = food_quality * min(e₀, (one(I) - σ) * iron_grazing_ratio)

    return total_specific_grazing * I, growth_efficiency
end

@inline function grazing(zoo::QualityDependantZooplankton, T, I, food_availability::NamedTuple, iron_availability::NamedTuple)
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensitivity
    p    = zoo.food_preferences
    food = keys(p)
    J    = zoo.specific_food_threshold_concentration
    K    = zoo.grazing_half_saturation
    food_threshold_concentration = zoo.food_threshold_concentration

    N = length(food)

    base_grazing_rate = g₀ * b ^ T

    total_food = sum(ntuple(n -> getproperty(food_availability, food[n]) * p[n], Val(N)))

    available_total_food = sum(ntuple(n -> max(zero(I), getproperty(food_availability, food[n]) - J) * p[n], Val(N)))

    concentration_limited_grazing = max(zero(I), available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food)

    θFe = zoo.iron_ratio
    e₀  = zoo.minimum_growth_efficiency
    σ   = zoo.non_assimilated_fraction

    total_iron = sum(ntuple(n -> getproperty(iron_availability, food[n]) * p[n], Val(N)))

    iron_grazing_ratio = total_iron / (θFe * total_specific_grazing + eps(zero(I)))

    food_quality = min(one(I), iron_grazing_ratio)

    growth_efficiency = food_quality * min(e₀, (one(I) - σ) * iron_grazing_ratio)

    return total_specific_grazing * I, growth_efficiency
end

@inline function grazing(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    # food quantity
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensitivity
    p    = zoo.food_preferences
    food = prey_names(bgc, val_name)
    J    = zoo.specific_food_threshold_concentration
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

    # food quality
    θFe = zoo.iron_ratio
    e₀  = zoo.minimum_growth_efficiency
    σ   = zoo.non_assimilated_fraction

    iron_availability = extract_iron_availability(bgc, i, j, k, fields, food)

    total_iron = sum(ntuple(n->iron_availability[n] * p[n], Val(N)))

    iron_grazing_ratio = total_iron / (θFe * total_specific_grazing + eps(0.0))

    food_quality = min(1, iron_grazing_ratio)

    growth_efficiency = food_quality * min(e₀, (1 - σ) * iron_grazing_ratio)

    return total_specific_grazing * I, growth_efficiency
end

@inline function flux_feeding(g₀,
                              b,
                              T,
                              I,
                              sinking_flux)

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_flux

    return total_specific_flux_feeding * I
end

@inline function flux_feeding(zoo::QualityDependantZooplankton, T, I, sinking_flux)
    return flux_feeding(zoo.maximum_flux_feeding_rate,
                        zoo.temperature_sensitivity,
                        T,
                        I,
                        sinking_flux)
end

@inline function flux_feeding(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = edible_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    return flux_feeding(zoo, T, I, sinking_flux)
end

@inline function mortality(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    b  = zoo.temperature_sensitivity
    m₀ = zoo.quadratic_mortality
    Kₘ = zoo.mortality_half_saturation
    r  = zoo.linear_mortality

    I = zooplankton_concentration(val_name, i, j, k, fields)
    T = @inbounds fields.T[i, j, k]

    O₂ = @inbounds fields.O₂[i, j, k]

    temperature_factor = b^T

    concentration_factor = I / (I + Kₘ)

    return temperature_factor * I * (m₀ * I + r * (concentration_factor + 3 * anoxia_factor(bgc, O₂)))
end

@inline function linear_mortality(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    b  = zoo.temperature_sensitivity
    Kₘ = zoo.mortality_half_saturation
    r  = zoo.linear_mortality
    
    T = @inbounds fields.T[i, j, k]
    O₂ = @inbounds fields.O₂[i, j, k]
    I = zooplankton_concentration(val_name, i, j, k, fields)

    temperature_factor = b^T

    concentration_factor = I / (I + Kₘ)

    return temperature_factor * r * (concentration_factor + 3 * anoxia_factor(bgc, O₂)) * I
end

#####
##### Effect on other compartments
#####

@inline function grazing(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensitivity
    p    = zoo.food_preferences
    food = prey_names(bgc, val_name)
    J    = zoo.specific_food_threshold_concentration
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

    P = concentration(val_prey_name, i, j, k, fields)

    return grazing_preference(val_prey_name, p) * max(0, P - J) * total_specific_grazing / (available_total_food + eps(0.0)) * I
end

@inline function flux_feeding(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀ =  zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensitivity

    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = flux_rate(val_prey_name, i, j, k, grid, fields, auxiliary_fields)

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_flux 

    return total_specific_flux_feeding * I
end

include("grazing_waste.jl")
include("mortality_waste.jl")
