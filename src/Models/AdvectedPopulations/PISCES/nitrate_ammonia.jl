"""
    NitrateAmmonia

A parameterisation for the evolution of nitrate (`NO₃`) and ammonia (`NH₄`)
where ammonia can be `nitrif`ied into nitrate, nitrate and ammonia are supplied
by the bacterial degredation of dissolved organic matter, and consumed by 
phytoplankton. Additionally waste produces ammonia through various means.

"""
@kwdef struct NitrateAmmonia{FT}
               maximum_nitrifcation_rate :: FT = 0.05 / day  # 1 / s
                   maximum_fixation_rate :: FT = 0.013 / day # mmol N / m³ (maybe shouldn't be a rate)
       iron_half_saturation_for_fixation :: FT = 0.1         # μmol Fe / m³
  phosphate_half_saturation_for_fixation :: FT = 0.8         # mmol P / m³
           light_saturation_for_fixation :: FT = 50.0        # W / m²
end

@inline function (nitrogen::NitrateAmmonia)(::Val{:NO₃}, bgc,
                                            x, y, z, t,
                                            P, D, Z, M, 
                                            PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, 
                                            SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, 
                                            CaCO₃, DIC, Alk, 
                                            O₂, T, S,
                                            zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(nitrogen, bgc, NH₄, O₂, mixed_layer_PAR) * θ

    remin = oxic_remineralisation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ) * θ

    nanophytoplankton_consumption = nitrate_uptake(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption = nitrate_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption) * θ

    return nitrif + remin - consumption # an extra term is present in Aumount 2015 but I suspect it is a typo
    # to conserve nitrogen I've dropped some ratios for denit etc, and now have bacterial_degregation go to denit in NO3 and remineralisation in NH4_half_saturation_const_for_DOC_remin
    # need to check...
end

@inline function (nitrogen::NitrateAmmonia)(::Val{:NH₄}, bgc,
                                            x, y, z, t,
                                            P, D, Z, M, 
                                            PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, 
                                            SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, 
                                            CaCO₃, DIC, Alk, 
                                            O₂, T, S,
                                            zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)
    θ = bgc.nitrogen_redfield_ratio

    nitrif = nitrification(nitrogen, bgc, NH₄, O₂, mixed_layer_PAR) * θ

    respiration_product = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T) * θ

    microzooplankton_grazing_waste = specific_inorganic_grazing_waste(bgc.microzooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, T, wPOC, wGOC) * Z
    mesozooplankton_grazing_waste  = specific_inorganic_grazing_waste(bgc.mesozooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, T, wPOC, wGOC) * M

    grazing_waste = (microzooplankton_grazing_waste + mesozooplankton_grazing_waste) * θ

    denit = denitrifcation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ) * θ

    nanophytoplankton_consumption = ammonia_uptake(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption = ammonia_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption) * θ

    fixation = nitrogen_fixation(nitrogen, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, PAR)

    # again an extra term is present in Aumount 2015 but I suspect it is a typo
    return fixation + respiration_product + grazing_waste + denit - consumption - nitrif
end

@inline function nitrification(nitrogen, bgc, NH₄, O₂, PAR)
    λ = nitrogen.maximum_nitrifcation_rate

    ΔO₂ = anoxia_factor(bgc, O₂)

    return λ * NH₄ / (1 + PAR) * (1 - ΔO₂)
end

@inline function nitrogen_fixation(nitrogen, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, PAR)
    Nₘ    = nitrogen.maximum_fixation_rate
    K_Fe  = nitrogen.iron_half_saturation_for_fixation
    K_PO₄ = nitrogen.phosphate_half_saturation_for_fixation
    E     = nitrogen.light_saturation_for_fixation

    phyto = bgc.nanophytoplankton

    _, _, _, LN, _, _ = phyto.nutrient_limitation(bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    fixation_limit = ifelse(LN >= 0.8, 0.01, 1 - LN)

    μ = base_production_rate(bgc.nanophytoplankton.growth_rate, T)

    growth_requirment = max(0, μ - 2.15)

    nutrient_limitation = min(Fe / (Fe +  K_Fe), PO₄ / (PO₄ + K_PO₄))

    light_limitation = 1 - exp(-PAR / E)

    return Nₘ * growth_requirment * fixation_limit * nutrient_limitation * light_limitation
end
