struct Phosphate end

@inline function (phosphate::Phosphate)(::Val{:PO₄}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi,
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, 
                                        zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    θ = bgc.phosphate_redfield_ratio

    microzooplankton_grazing_waste = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * Z
    mesozooplankton_grazing_waste  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * M

    grazing_waste = microzooplankton_grazing_waste + mesozooplankton_grazing_waste

    respiration_product = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    remineralisation = oxic_remineralisation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    denit = denitrifcation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    nanophytoplankton_consumption = total_growth(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption = total_growth(bgc.diatoms, D, DChl, DsFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption)

    return θ * (grazing_waste + respiration_product + remineralisation + denit - consumption)
end
