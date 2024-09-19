"""
    CarbonateSystem

Default parameterisation for `DIC`` and `Alk`alinity evolution. 
"""
struct CarbonateSystem end

@inline function (carbonates::CarbonateSystem)(::Val{:DIC}, bgc,
                                               x, y, z, t,
                                               P, D, Z, M, 
                                               PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, 
                                               SFe, BFe, PSi, 
                                               NO₃, NH₄, PO₄, Fe, Si, 
                                               CaCO₃, DIC, Alk, 
                                               O₂, T, S,
                                               zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    microzooplankton_respiration = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * Z
    mesozooplankton_respiration  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * M

    zooplankton_respiration = microzooplankton_respiration + mesozooplankton_respiration

    upper_trophic_respiration = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    dissolved_degredation = bacterial_degradation(bgc.dissolved_organic_matter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    calcite_diss = calcite_dissolution(bgc.calcite, CaCO₃, Ω)

    calcite_prod = calcite_production(bgc.calcite, bgc, z, P, D, PChl, PFe, Z, M, POC, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_consumption, = total_production(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    diatom_consumption,            = total_production(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = nanophytoplankton_consumption + diatom_consumption

    return zooplankton_respiration + upper_trophic_respiration + dissolved_degredation + calcite_diss - calcite_prod - consumption
end

@inline function (carbonates::CarbonateSystem)(::Val{:Alk}, bgc, args...)
    θ = bgc.nitrogen_redfield_ratio

    nitrate_production = bgc.nitrogen(Val(:NO₃), bgc, args...) * θ
    ammonia_production = bgc.nitrogen(Val(:NH₄), bgc, args...) * θ
    calcite_production = bgc.calcite(Val(:CaCO₃), bgc, args...)

    # I think there are typos in Aumount 2015 but this is what it should be

    return ammonia_production - nitrate_production - 2 * calcite_production
end
