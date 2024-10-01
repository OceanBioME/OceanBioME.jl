"""
    NitrateAmmonia

A parameterisation for the evolution of nitrate (`NO₃`) and ammonia (`NH₄`)
where ammonia can be `nitrif`ied into nitrate, nitrate and ammonia are supplied
by the bacterial degredation of dissolved organic matter, and consumed by 
phytoplankton. Additionally waste produces ammonia through various means.

"""
@kwdef struct NitrateAmmonia{FT}
              maximum_nitrification_rate :: FT = 0.05  / day # 1 / s
                   maximum_fixation_rate :: FT = 0.013 / day # mmol N / m³ (maybe shouldn't be a rate)
       iron_half_saturation_for_fixation :: FT = 0.1         # μmol Fe / m³
  phosphate_half_saturation_for_fixation :: FT = 0.8         # mmol P / m³
           light_saturation_for_fixation :: FT = 50.0        # W / m²
end

required_biogeochemical_tracers(::NitrateAmmonia) = (:NO₃, :NH₄)

const NitrateAmnmoniaPISCES = PISCES{<:Any, <:Any, <:Any, <:Any, <:NitrateAmmonia}

@inline function (bgc::NitrateAmnmoniaPISCES)(i, j, k, grid, val_name::Val{:NO₃}, clock, fields, auxiliary_fields)
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(bgc.nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    remineralisation = oxic_remineralisation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    consumption = uptake(bgc.phytoplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return nitrif + θ * (remineralisation - consumption)
end

@inline function (bgc::NitrateAmnmoniaPISCES)(i, j, k, grid, val_name::Val{:NH₄}, clock, fields, auxiliary_fields)
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(bgc.nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    remineralisation = anoxic_remineralisation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    consumption = uptake(bgc.phytoplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    grazing_waste = inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic_waste = upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    fixation = nitrogen_fixation(bgc.nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return fixation + θ * (remineralisation + grazing_waste + upper_trophic_waste - consumption) - nitrif
end

@inline function nitrification(nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    λ = nitrogen.maximum_nitrification_rate

    O₂  = @inbounds fields.O₂[i, j, k]
    NH₄ = @inbounds fields.NH₄[i, j, k]

    PAR = @inbounds auxiliary_fields.mixed_layer_PAR[i, j, k]

    ΔO₂ = anoxia_factor(bgc, O₂)

    return λ * NH₄ / (1 + PAR) * (1 - ΔO₂)
end

@inline function nitrogen_fixation(nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    Nₘ    = nitrogen.maximum_fixation_rate
    K_Fe  = nitrogen.iron_half_saturation_for_fixation
    K_PO₄ = nitrogen.phosphate_half_saturation_for_fixation
    E     = nitrogen.light_saturation_for_fixation

    PAR = @inbounds auxiliary_fields.PAR[i, j, k]

    Fe  = @inbounds  fields.Fe[i, j, k]
    PO₄ = @inbounds fields.PO₄[i, j, k]

    availability_limitation = nitrogen_availability_limitation(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    fixation_limit = ifelse(availability_limitation >= 0.8, 0.01, 1 - availability_limitation)

    μ = base_production_rate(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    growth_requirment = max(0, μ - 2.15)

    nutrient_limitation = min(Fe / (Fe +  K_Fe), PO₄ / (PO₄ + K_PO₄))

    light_limitation = 1 - exp(-PAR / E)

    return Nₘ * growth_requirment * fixation_limit * nutrient_limitation * light_limitation
end