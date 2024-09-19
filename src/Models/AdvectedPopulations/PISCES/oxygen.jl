@kwdef struct Oxygen{FT}
     ratio_for_respiration :: FT = 133/122 # mol O₂ / mol C
    ratio_for_nitrifcation :: FT = 32/122  # mol O₂ / mol C
end

@inline function (oxy::Oxygen)(::Val{:O₂}, bgc,
                               x, y, z, t,
                               P, D, Z, M, 
                               PChl, DChl, PFe, DFe, DSi,
                               DOC, POC, GOC, 
                               SFe, BFe, PSi, 
                               NO₃, NH₄, PO₄, Fe, Si, 
                               CaCO₃, DIC, Alk, 
                               O₂, T, S,
                               zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    θ_resp = oxy.ratio_for_respiration
    θ_nitrif  = oxy.ratio_for_nitrifcation

    microzooplankton_respiration = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * Z
    mesozooplankton_respiration  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * M

    zooplankton_respiration = θ_resp * (microzooplankton_respiration + mesozooplankton_respiration)

    upper_trophic_respiration = θ_resp * inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    remin = (θ_resp + θ_nitrif) * oxic_remineralisation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ)
    denit = θ_resp * denitrifcation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ)

    nitrif = θ_nitrif * nitrification(bgc.nitrogen, bgc, NH₄, O₂, mixed_layer_PAR) 
    
    nanophytoplankton_nitrate_consumption = nitrate_uptake(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_nitrate_consumption = nitrate_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    nitrate_consumption = (θ_resp + θ_nitrif) * (nanophytoplankton_nitrate_consumption + diatom_nitrate_consumption)

    nanophytoplankton_ammonia_consumption = ammonia_uptake(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_ammonia_consumption = ammonia_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    ammonia_consumption = θ_resp * (nanophytoplankton_ammonia_consumption + diatom_ammonia_consumption)

    photosynthesis = nitrate_consumption + ammonia_consumption

    fixation = θ_nitrif * nitrogen_fixation(bgc.nitrogen, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, PAR)

    return photosynthesis + fixation - zooplankton_respiration - upper_trophic_respiration - nitrif - remin - denit
end
