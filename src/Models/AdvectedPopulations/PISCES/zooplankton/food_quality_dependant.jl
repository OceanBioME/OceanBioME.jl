"""
    QualityDependantZooplankton

The PISCES zooplankton growth model where each class has preferences
for grazing on nanophytoplankton (P), diatoms (D), microzooplankton (Z),
and particulate organic matter (POC), and can flux feed on sinking 
particulates (POC and GOC).

This model assumes a fixed ratio for all other elements (i.e. N, P, Fe).
"""
@kwdef struct QualityDependantZooplankton{FT, FP}
                   temperature_sensetivity :: FT = 1.079 #
                      maximum_grazing_rate :: FT         # 1 / s

                          food_preferences :: FP 

              food_threshold_concentration :: FT = 0.3   # mmol C / m³
    specific_food_thresehold_concentration :: FT = 0.001 # mmol C / m³

                   grazing_half_saturation :: FT = 20.0  # mmol C / m³

                 maximum_flux_feeding_rate :: FT         # m / (mmol C / m³)

                                iron_ratio :: FT         # μmol Fe / mmol C

                 minimum_growth_efficiency :: FT         #
                  non_assililated_fraction :: FT = 0.3   #

                 mortality_half_saturation :: FT = 0.2   # mmol C / m³
                       quadratic_mortality :: FT         # 1 / (mmol C / m³) / s
                          linear_mortality :: FT         # 1 / s

    # this should be called inorganic excretion factor
              dissolved_excretion_fraction :: FT = 0.6   #
              undissolved_calcite_fraction :: FT         #
end

required_biogeochemical_tracers(::QualityDependantZooplankton, name_base) = tuple(name_base)

@inline function growth_death(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    gI, e = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 
    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 
    mI = mortality(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return e * (gI + gfI) - mI
end

@inline extract_food_availability(i, j, k, fields, names::NTuple{N}) where N =
    ntuple(n -> concentration(Val(names[n]), i, j, k, fields), Val(N))

@inline extract_iron_availability(i, j, k, bgc, fields, names::NTuple{N}) where N =
    ntuple(n -> iron_ratio(Val(names[n]), i, j, k, bgc, fields), Val(N))

@inline function grazing(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
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

    food_availability = extract_food_availability(i, j, k, fields, food)

    total_food = sum(ntuple(n->food_availability[n] * p[n], Val(N)))

    available_total_food = sum(ntuple(n->max(zero(grid), (food_availability[n] - J)) * p[n], Val(N)))

    concentration_limited_grazing = max(0, available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food)

    # food quality
    θFe = zoo.iron_ratio
    e₀  = zoo.minimum_growth_efficiency
    σ   = zoo.non_assililated_fraction

    iron_availabillity = extract_iron_availability(i, j, k, bgc, fields, food)

    total_iron = sum(ntuple(n->iron_availabillity[n] * p[n], Val(N)))

    iron_grazing_ratio = total_iron / (θFe * total_specific_grazing + eps(0.0))

    food_quality = min(1, iron_grazing_ratio)

    growth_efficiency = food_quality * min(e₀, (1 - σ) * iron_grazing_ratio)

    return total_specific_grazing * I, growth_efficiency
end

@inline function flux_feeding(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀ =  zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensetivity

    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = edible_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_flux 

    return total_specific_flux_feeding * I
end

@inline function mortality(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    b  = zoo.temperature_sensetivity
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
    b  = zoo.temperature_sensetivity
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
##### Effect on other compartements
#####

@inline function grazing(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensetivity
    p    = zoo.food_preferences
    food = zoo.food_names
    J    = zoo.specific_food_thresehold_concentration
    K    = zoo.grazing_half_saturation
    food_threshold_concentration = zoo.food_threshold_concentration

    N = length(food)

    I = zooplankton_concentration(val_name, i, j, k, fields)
    T = @inbounds fields.T[i, j, k]

    base_grazing_rate = g₀ * b ^ T

    food_availability = extract_food_availability(i, j, k, fields, food)

    total_food = sum(ntuple(n->food_availability[n] * p[n], Val(N)))

    available_total_food = sum(ntuple(n->max(zero(grid), (food_availability[n] - J)) * p[n], Val(N)))

    concentration_limited_grazing = max(0, available_total_food - min(available_total_food / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + total_food) 

    P = concentration(val_prey_name, i, j, k, fields)

    return grazing_preference(val_prey_name, p) * max(0, P - J) * total_specific_grazing / (available_total_food + eps(0.0)) * I
end

@inline function flux_feeding(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    g₀ =  zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensetivity

    I = zooplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    sinking_flux = flux_rate(val_prey_name, i, j, k, grid, fields, auxiliary_fields)

    base_flux_feeding_rate = g₀ * b ^ T

    total_specific_flux_feeding = base_flux_feeding_rate * sinking_flux 

    return total_specific_flux_feeding * I
end

include("grazing_waste.jl")
include("mortality_waste.jl")