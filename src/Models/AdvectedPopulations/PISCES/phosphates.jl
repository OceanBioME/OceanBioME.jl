"""
    Phosphate

Evolution of phosphate (PO₄).
"""
struct Phosphate end

@inline function (phosphate::Phosphate)(::Val{:PO₄}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi,
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

    θ = bgc.phosphate_redfield_ratio

    microzooplankton_grazing_waste = specific_inorganic_grazing_waste(bgc.microzooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, T, wPOC, wGOC) * Z
    mesozooplankton_grazing_waste  = specific_inorganic_grazing_waste(bgc.mesozooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, T, wPOC, wGOC) * M

    grazing_waste = microzooplankton_grazing_waste + mesozooplankton_grazing_waste

    respiration_product = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    remineralisation = oxic_remineralisation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ)

    denit = denitrifcation(bgc.dissolved_organic_matter, bgc, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, O₂, T, zₘₓₗ, zₑᵤ)

    nanophytoplankton_consumption, = total_production(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    diatom_consumption, = total_production(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = (nanophytoplankton_consumption + diatom_consumption)

    return θ * (grazing_waste + respiration_product + remineralisation + denit - consumption)
end
