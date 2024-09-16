@kwdef struct TwoCompartementParticulateOrganicMatter{FT, AP}
                      temperature_sensetivity :: FT = 1.066
                          base_breakdown_rate :: FT = 0.025 / day
                       aggregation_parameters :: AP = (0, 0, 0, 0)#(25.9, 4452, 3.3, 47.1) .* (10^-6 / day)
                 minimum_iron_scavenging_rate :: FT = 0.0#3e-5/day
           load_specific_iron_scavenging_rate :: FT = 0.0#0.005/day
           dust_specific_iron_scavenging_rate :: FT = 0.0#150/day
  small_fraction_of_bacterially_consumed_iron :: FT = 0.5
  large_fraction_of_bacterially_consumed_iron :: FT = 0.5 
                base_liable_silicate_fraction :: FT = 0.5
            fast_dissolution_rate_of_silicate :: FT = 0.025/day
            slow_dissolution_rate_of_silicate :: FT = 0.003/day
end

@inline function specific_degredation_rate(poc::TwoCompartementParticulateOrganicMatter, bgc, O₂, T)
    λ₀ = poc.base_breakdown_rate
    b  = poc.temperature_sensetivity

    ΔO₂ = anoxia_factor(bgc, O₂)

    return λ₀ * b^T * (1 - 0.45 * ΔO₂)
end

@inline function aggregation(poc::TwoCompartementParticulateOrganicMatter, bgc, z, POC, GOC, zₘₓₗ)
    a₁, a₂, a₃, a₄ = poc.aggregation_parameters

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    return shear * (a₁ * POC^2 + a₂ * POC * GOC) + a₃ * POC * GOC + a₄ * POC^2
end

include("particulate_organic_carbon.jl")
include("iron_in_particles.jl")
include("silicon_in_particles.jl")