@kwdef struct Oxygen{FT}
     ratio_for_respiration :: FT = 133/122 # I think this is a more accurate name
    ratio_for_nitrifcation :: FT = 32/122
end

@inline function (oxy::Oxygen)(::Val{:O₂}, bgc,
                               x, y, z, t,
                               P, D, Z, M, 
                               PChl, DChl, PFe, DFe, DSi,
                               DOC, POC, GOC, 
                               SFe, BFe, PSi, 
                               NO₃, NH₄, PO₄, Fe, Si, 
                               CaCO₃, DIC, Alk, 
                               O₂, T, 
                               zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    θ_resp = oxy.ratio_for_respiration
    θ_nitrif  = oxy.ratio_for_nitrifcation

    microzooplankton_respiration = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * Z
    mesozooplankton_respiration  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * M

    zooplankton_respiration = θ_resp * (microzooplankton_respiration + mesozooplankton_respiration)

    upper_trophic_respiration = θ_resp * inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    # not sure I've got this right
    remin = θ_resp * oxic_remineralisation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    nitrif = θ_nitrif * nitrification(nitrogen, NH₄, O₂, mixed_layer_PAR) 
    
    nanophytoplankton_nitrate_consumption = nitrate_uptake(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_nitrate_consumption = nitrate_uptake(bgc.diatoms, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    nitrate_uptake = (θ_resp + θ_nitrif) * (nanophytoplankton_nitrate_consumption + diatom_nitrate_consumption)

    nanophytoplankton_ammonia_consumption = ammonia_uptake(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_ammonia_consumption = ammonia_uptake(bgc.diatoms, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    ammonia_uptake = θ_resp * (nanophytoplankton_ammonia_consumption + diatom_ammonia_consumption)

    photosynthesis = nitrate_uptake + ammonia_uptake

    return photosynthesis - zooplankton_respiration - upper_trophic_respiration - nitrif - remin
end
