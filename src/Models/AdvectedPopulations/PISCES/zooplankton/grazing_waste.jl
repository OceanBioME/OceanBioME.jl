include("iron_grazing.jl")

@inline function non_assimilated_waste(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    σ = zoo.non_assililated_fraction

    gI, = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    return σ * (gI + gfI)
end

@inline function excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    σ = zoo.non_assililated_fraction

    gI, growth_efficiency = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    return (1 - σ - e) * (gI + gfI)
end

@inline function inorganic_excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    γ = zoo.dissolved_excretion_fraction

    return γ * excretion(zoo, val_name, i, j, k, grid, bgc, clock, fields)
end

@inline function organic_excretion(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    γ = zoo.dissolved_excretion_fraction

    return (1 - γ) * excretion(zoo, val_name, i, j, k, grid, bgc, clock, fields)
end

@inline function non_assimilated_iron_waste(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    σ = zoo.non_assililated_fraction

    gI = iron_grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    gfI = iron_flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    return σ * (gI + gfI)
end

@inline function non_assimilated_iron(zoo::QualityDependantZooplankton, val_name, i, j, k, grid, bgc, clock, fields)
    θ = zoo.iron_ratio
    σ = zoo.non_assililated_fraction

    gI, growth_efficiency = grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields) 
    gfI = flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields) 

    total_carbon = gI + gfI

    total_iron = (iron_grazing(zoo, val_name, i, j, k, grid, bgc, clock, fields) 
                  + iron_flux_feeding(zoo, val_name, i, j, k, grid, bgc, clock, fields))

    grazing_iron_ratio = (1 - σ) * total_iron / (total_carbon + eps(0.0))

    non_assimilated_iron_ratio = max(0, grazing_iron_ratio - growth_efficiency * θ)

    return non_assililated_fraction * gI
end

@inline function calcite_loss(zoo::QualityDependantZooplankton, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields)
    η = zoo.undissolved_calcite_fraction
    
    g = grazing(zoo, val_name, val_prey_name, i, j, k, grid, bgc, clock, fields)

    return η * g
end
