@kwdef struct NitrateAmmonia{FT}
               maximum_nitrifcation_rate :: FT = 0.05 / day
    anoxic_denitrifcation_nitrogen_ratio :: FT = 0.86
                   maximum_fixation_rate :: FT = 0.013 / day
       iron_half_saturation_for_fixation :: FT = 0.1
  phosphate_half_saturation_for_fixation :: FT = 0.8
           light_saturation_for_fixation :: FT = 50.0
end

@inline function (nitrogen::NitrateAmmonia)(bgc, val_name::Val{:NO₃}, 
                                            x, y, z, t,
                                            P, D, Z, M, 
                                            PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, 
                                            SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, 
                                            CaCO₃, DIC, Alk, 
                                            O₂, T, 
                                            zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)
    R = nitrogen.anoxic_denitrifcation_nitrogen_ratio
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(nitrogen, NH₄, O₂, mixed_layer_PAR) * θ

    denit = R * denitrifcation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    nanophytoplankton_consumption = nitrate_uptake(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption = nitrate_uptake(bgc.diatoms, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption) * θ

    return nitrif - denit - consumption # an extra term is present in Aumount 2015 but I suspect it is a typo
end

@inline function (nitrogen::NitrateAmmonia)(bgc, val_name::Val{:NH₄}, 
                                            x, y, z, t,
                                            P, D, Z, M, 
                                            PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, 
                                            SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, 
                                            CaCO₃, DIC, Alk, 
                                            O₂, T, 
                                            zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)
    R = nitrogen.anoxic_denitrifcation_nitrogen_ratio
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(nitrogen, NH₄, O₂, mixed_layer_PAR) * θ

    denit = R * denitrifcation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    respiration_product = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T) * θ

    microzooplankton_grazing_waste = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * Z
    mesozooplankton_grazing_waste  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * M

    grazing_waste = (microzooplankton_grazing_waste + mesozooplankton_grazing_waste) * θ

    remineralisation = oxic_remineralisation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ) * θ

    nanophytoplankton_consumption = ammonia_uptake(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption = ammonia_uptake(bgc.diatoms, D, DChl, DsFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption) * θ

    fixation = nitrogen_fixation(nitrogen, bgc, NO₃, NH₄, PO₄, Fe, Si, Si′, PAR)

    # again an extra term is present in Aumount 2015 but I suspect it is a typo
    return fixation + respiration_product + grazing_waste + denit + remineralisation - consumption - nitrif
end

@inline function nitrification(nitrogen, NH₄, O₂, PAR)
    λ = nitrogen.maximum_nitrifcation_rate

    ΔO₂ = anoxia_factor(bgc, O₂)

    return λ * NH₄ / (1 + PAR) * (1 - ΔO₂)
end

@inline function nitrogen_fixation(nitrogen, bgc, NO₃, NH₄, PO₄, Fe, Si, Si′, PAR)
    Nₘ    = nitogen.maximum_fixation_rate
    K_Fe  = nitrogen.iron_half_saturation_for_fixation
    K_PO₄ = nitrogen.phosphate_half_saturation_for_fixation
    E     = nitrogen.light_saturation_for_fixation

    phyto = bgc.nanophytoplankton

    _, _, _, LN, _, _ = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    fixation_limit = ifelse(LN >= 0.8, 0.01, 1 - LN)

    μ = base_production_rate(bgc.nanophytoplankton.growth_rate, T)

    growth_requirment = max(0, μ - 2.15)

    nutrient_limitation = min(Fe / (Fe +  K_Fe), PO₄ / (PO₄ + K_PO₄))

    light_limitation = 1 - exp(-PAR / E)

    return Nₘ * growth_requirment * fixation_limit * nutrient_limitation * light_limitation
end
