struct CarbonateSystem end

@inline function (carbonates::CarbonateSystem)(bgc, ::Val{:DIC}, 
                                               x, y, z, t,
                                               P, D, Z, M, 
                                               PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, 
                                               SFe, BFe, PSi, 
                                               NO₃, NH₄, PO₄, Fe, Si, 
                                               CaCO₃, DIC, Alk, 
                                               O₂, T, 
                                               zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    microzooplankton_respiration = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * Z
    mesozooplankton_respiration  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, P, D, PFe, DFe, Z, POC, SFe) * M

    zooplankton_respiration = microzooplankton_respiration + mesozooplankton_respiration

    upper_trophic_respiration = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    dissolved_degredation = bacterial_degradation(bgc.dissolved_organic_matter, z, Z, M, DOM, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    calcite_diss = calcite_dissolution(bgc.calcite, CaCO₃, Ω)

    calcite_prod = calcite_production(bgc.calcite, bgc, P, D, PChl, PFe, Z, M, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    calcite = calcite_diss - calcite_prod

    nanophytoplankton_consumption = total_production(bgc.nanophytoplankton, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    diatom_consumption            = total_production(bgc.diatoms, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = nanophytoplankton_consumption + diatom_consumption

    return zooplankton_respiration + upper_trophic_respiration + dissolved_degredation + calcite - consumption
end

@inline function (carbonates::CarbonateSystem)(bgc, ::Val{:Alk}, args...)
    nitrate_production = bgc.nitrogen(bgc, Val(:NO₃), args...)
    ammonia_production = bgc.nitrogen(bgc, Val(:NO₃), args...)
    calcite_production = bgc.calcite(bgc, Val(:CaCO₃), args...)

    # I think there are typos in Aumount 2015 but this is what it should be
    return ammonia_production - nitrate_production - 2 * calcite_production
end
