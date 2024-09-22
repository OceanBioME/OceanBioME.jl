"""
    TwoCompartementParticulateOrganicMatter

A quota parameterisation for particulate organic matter with two size classes,
each with carbon and iron compartements, and a silicate compartement for the
large size class.

Confusingly we decided to name these compartmenets `POC` and `GOC` for the small
and large carbon classes, `SFe` and `BFe` for the small and  ̶l̶a̶r̶g̶e̶ big iron 
compartements, and `PSi` for the  ̶l̶a̶r̶g̶e̶ particulate silicon (*not* the 
phytoplankton silicon).
"""
@kwdef struct TwoCompartementParticulateOrganicMatter{FT, AP}
                      temperature_sensetivity :: FT = 1.066        #
                          base_breakdown_rate :: FT = 0.025 / day  # 1 / s
# (1 / (mmol C / m³),  1 / (mmol C / m³),  1 / (mmol C / m³) / s,  1 / (mmol C / m³) / s)
                       aggregation_parameters :: AP = (25.9, 4452, 3.3, 47.1) .* (10^-6 / day)
                 minimum_iron_scavenging_rate :: FT = 3e-5/day     # 1 / s
           load_specific_iron_scavenging_rate :: FT = 0.005/day    # 1 / (mmol C / m³) / s
  small_fraction_of_bacterially_consumed_iron :: FT = 0.5          #
  large_fraction_of_bacterially_consumed_iron :: FT = 0.5          #
                base_liable_silicate_fraction :: FT = 0.5          #
            fast_dissolution_rate_of_silicate :: FT = 0.025/day    # 1 / s
            slow_dissolution_rate_of_silicate :: FT = 0.003/day    # 1 / s
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