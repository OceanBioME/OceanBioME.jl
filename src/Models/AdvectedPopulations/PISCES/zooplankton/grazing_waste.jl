include("iron_grazing.jl")

@inline function non_assimilated_waste(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    σ = zoo.non_assimilated_fraction

    gI, = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return σ * (gI + gfI)
end

@inline function excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    σ = zoo.non_assimilated_fraction

    gI, e = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return (1 - σ - e) * (gI + gfI)
end

@inline function inorganic_excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    γ = zoo.dissolved_excretion_fraction

    return γ * excretion(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
end

@inline function organic_excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    γ = zoo.dissolved_excretion_fraction

    return (1 - γ) * excretion(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
end

@inline function non_assimilated_iron_waste(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    σ = zoo.non_assimilated_fraction

    gI = iron_grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    gfI = iron_flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return σ * (gI + gfI)
end

@inline function non_assimilated_iron(θ,
                                      σ,
                                      g₀,
                                      b,
                                      pP,
                                      pD,
                                      pZ,
                                      pPOC,
                                      J,
                                      K,
                                      food_threshold_concentration,
                                      e₀,
                                      gᶠ₀,
                                      T,
                                      I,
                                      food_availability::NamedTuple,
                                      iron_availability::NamedTuple,
                                      sinking_flux,
                                      sinking_iron_flux)

    gI, growth_efficiency = grazing(g₀,
                                    b,
                                    pP,
                                    pD,
                                    pZ,
                                    pPOC,
                                    J,
                                    K,
                                    food_threshold_concentration,
                                    θ,
                                    e₀,
                                    σ,
                                    T,
                                    I,
                                    food_availability,
                                    iron_availability)
    gfI = flux_feeding(gᶠ₀,
                       b,
                       T,
                       I,
                       sinking_flux)

    zoo_assimilated_iron = θ * growth_efficiency * (gI + gfI)

    gIFe = iron_grazing(g₀,
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
                        food_availability,
                        iron_availability)

    gfIFe = iron_flux_feeding(gᶠ₀,
                              b,
                              T,
                              I,
                              sinking_iron_flux)

    lost_to_particles = σ * (gIFe + gfIFe)

    total_iron_grazed = gIFe + gfIFe

    return total_iron_grazed - lost_to_particles - zoo_assimilated_iron
end

@inline function non_assimilated_iron(zoo::QualityDependantZooplankton,
                                         T,
                                         I,
                                         food_availability::NamedTuple,
                                         iron_availability::NamedTuple,
                                         sinking_flux,
                                         sinking_iron_flux)
    return non_assimilated_iron(zoo.iron_ratio,
                                zoo.non_assimilated_fraction,
                                zoo.maximum_grazing_rate,
                                zoo.temperature_sensitivity,
                                zoo.food_preferences.P,
                                zoo.food_preferences.D,
                                zoo.food_preferences.Z,
                                zoo.food_preferences.POC,
                                zoo.specific_food_threshold_concentration,
                                zoo.grazing_half_saturation,
                                zoo.food_threshold_concentration,
                                zoo.minimum_growth_efficiency,
                                zoo.maximum_flux_feeding_rate,
                                T,
                                I,
                                food_availability,
                                iron_availability,
                                sinking_flux,
                                sinking_iron_flux)
end

@inline function non_assimilated_iron(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    θ = zoo.iron_ratio
    σ = zoo.non_assimilated_fraction

    gI, growth_efficiency = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 
    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    zoo_assimilated_iron = θ * growth_efficiency * (gI + gfI)

    gIFe = iron_grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    gfIFe = iron_flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    lost_to_particles = σ * (gIFe + gfIFe)

    total_iron_grazed = gIFe + gfIFe

    return total_iron_grazed - lost_to_particles - zoo_assimilated_iron # feels like a more straight forward way to write it
end

@inline function calcite_loss(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    η = zoo.undissolved_calcite_fraction
    
    g = grazing(zoo, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return η * g
end
