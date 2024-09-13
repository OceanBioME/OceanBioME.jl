@kwdef struct DissolvedOrganicMatter{FT, AP}
                             remineralisation_rate :: FT = 0.3/day
           microzooplankton_bacteria_concentration :: FT = 0.7
            mesozooplankton_bacteria_concentration :: FT = 1.4
                    maximum_bacteria_concentration :: FT = 4.0
             bacteria_concentration_depth_exponent :: FT = 0.684
                  reference_bacteria_concentration :: FT = 1.0
                           temperature_sensetivity :: FT = 1.066
        doc_half_saturation_for_bacterial_activity :: FT = 417.0
    nitrate_half_saturation_for_bacterial_activity :: FT = 0.03
    ammonia_half_saturation_for_bacterial_activity :: FT = 0.003
  phosphate_half_saturation_for_bacterial_activity :: FT = 0.003
       iron_half_saturation_for_bacterial_activity :: FT = 0.01
                            aggregation_parameters :: AP = (0.37, 102, 3530, 5095, 114) .* (10^-6 / day) #(μmolCL⁻¹)⁻¹s⁻¹
                    maximum_iron_ratio_in_bacteria :: FT = 10^-3
                 iron_half_saturation_for_bacteria :: FT = 0.03
                     maximum_bacterial_growth_rate :: FT = 0.6 / day
end

@inline function (dom::DissolvedOrganicMatter)(bgc, ::Val{:DOC}, 
                                               x, y, z, t,
                                               P, D, Z, M, 
                                               PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, 
                                               SFe, BFe, PSi, 
                                               NO₃, NH₄, PO₄, Fe, Si, 
                                               CaCO₃, DIC, Alk, 
                                               O₂, T, 
                                               zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    nanophytoplankton_exudation = dissolved_exudate(bgc.nanophytoplankton, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    diatom_exudation = dissolved_exudate(bgc.nanophytoplankton, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    phytoplankton_exudation = nanophytoplankton_exudation + diatom_exudation

    particulate_degredation = dissolved_degredation_product(bgc.particulate_organic_matter, POC, GOC, O₂, T)

    respiration_product = dissolved_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    microzooplankton_grazing_waste = specific_dissolved_grazing_waste(bgc.microzooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * Z
    mesozooplankton_grazing_waste  = specific_dissolved_grazing_waste(bgc.mesozooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * M

    grazing_waste = microzooplankton_grazing_waste + mesozooplankton_grazing_waste

    degredation = bacterial_degradation(dom, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    aggregation_to_particles = aggregation(dom, bgc, z, DOC, POC, GOC, zₘₓₗ)

    return phytoplankton_exudation + particulate_degredation + respiration_product + grazing_waste - degredation - aggregation_to_particles
end

@inline function bacteria_concentration(dom::DissolvedOrganicMatter, z, Z, M, zₘₓₗ, zₑᵤ)
    bZ = dom.microzooplankton_bacteria_concentration
    bM = dom.mesozooplankton_bacteria_concentration
    a  = dom.bacteria_concentration_depth_exponent

    zₘ = min(zₘₓₗ, zₑᵤ)

    surface_bacteria = min(4, bZ * Z + bM * M)

    depth_factor = (zₘ / z) ^ a

    return ifelse(z >= zₘ, 1, depth_factor) * surface_bacteria
end

@inline function bacteria_activity(dom::DissolvedOrganicMatter, DOC, NO₃, NH₄, PO₄, Fe)
    K_DOC = dom.doc_half_saturation_for_bacterial_activity
    K_NO₃ = dom.nitrate_half_saturation_for_bacterial_activity
    K_NH₄ = ammonia_half_saturation_for_bacterial_activity
    K_PO₄ = phosphate_half_saturation_for_bacterial_activity
    K_Fe  = iron_half_saturation_for_bacterial_activity

    DOC_limit = DOC / (DOC + K_DOC)

    L_N   = (K_NO₃ * NH₄ + K_NH₄ * NO₃) / (K_NO₃ * K_NH₄ + K_NO₃ * NH₄ + K_NH₄ * NO₃)

    L_PO₄ = PO₄ / (PO₄ + K_PO₄)

    L_Fe  = Fe / (Fe + K_Fe)

    # assuming typo in paper otherwise it doesn't make sense to formulate L_NH₄ like this
    limiting_quota = min(L_N, L_PO₄, L_Fe)

    return limiting_quota * DOC_limit
end

@inline function bacterial_degradation(dom::DissolvedOrganicMatter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)
    Bact_ref = dom.reference_bacteria_concentration
    b = dom.temperature_sensetivity
    λ = dom.remineralisation_rate

    f = b^T

    Bact = bacteria_concentration(dom, z, Z, M, zₘₓₗ, zₑᵤ)

    LBact = bacteria_activity(dom, DOC, NO₃, NH₄, PO₄, Fe)

    return λ * f * LBact * Bact / Bact_ref * DOM # differes from Aumont 2015 since the dimensions don't make sense 
end

@inline function oxic_remineralisation(dom::DissolvedOrganicMatter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)
    ΔO₂ = anoxia_factor(bgc, O₂)

    degredation = bacterial_degradation(dom, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    return (1 - ΔO₂) * degredation
end

@inline function denitrifcation(dom::DissolvedOrganicMatter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)
    ΔO₂ = anoxia_factor(bgc, O₂)

    degredation = bacterial_degradation(dom, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    return ΔO₂ * degredation
end

@inline function aggregation(dom::DissolvedOrganicMatter, bgc, z, DOC, POC, GOC, zₘₓₗ)
    a₁, a₂, a₃, a₄, a₅ = dom.aggregation_parameters

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    Φ₁ = shear * (a₁ * DOC + a₂ * POC) * DOC
    Φ₂ = shear * (a₃ * GOC) * DOC
    Φ₃ = (a₄ * POC + a₅ * DOC) * DOC

    return Φ₁ + Φ₂ + Φ₃, Φ₁, Φ₂, Φ₃
end

@inline function bacterial_iron_uptake(dom::DissolvedOrganicMatter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)
    μ₀ = dom.maximum_bacterial_growth_rate
    b  = dom.temperature_sensetivity
    θ  = dom.iron_half_saturation_for_bacteria
    K  = dom.iron_half_saturation_for_bacteria

    μ = μ₀ * b^T

    Bact = bacteria_concentration(dom, z, Z, M, zₘₓₗ, zₑᵤ)

    L = bacteria_activity(dom, DOC, NO₃, NH₄, PO₄, Fe)

    return μ * L * θ * Fe / (Fe + K) * Bact
end