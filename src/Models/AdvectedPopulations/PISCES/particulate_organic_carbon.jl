@inline function (poc::TwoCompartementParticulateOrganicMatter)(bgc, ::Val{:POC}, 
                                                                x, y, z, t,
                                                                P, D, Z, M, 
                                                                PChl, DChl, PFe, DFe, DSi,
                                                                DOC, POC, GOC, 
                                                                SFe, BFe, PSi, 
                                                                NO₃, NH₄, PO₄, Fe, Si, 
                                                                CaCO₃, DIC, Alk, 
                                                                O₂, T, 
                                                                zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    grazing_waste = specific_non_assimilated_waste(bgc.microzooplankton, bgc, P, D, Z, POC, GOC) * Z

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, D, zₘₓₗ)

    nanophytoplankton_mortality = (1 - 0.5 * R_CaCO₃) * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality)

    diatom_linear_mortality, = mortality(bgc.diatoms, bgc, z, D, zₘₓₗ)

    diatom_mortality = 0.5 * diatom_linear_mortality

    microzooplankton_mortality = mortality(bgc.microzooplankton, bgc, I, O₂, T)

    # degredation
    λ = specific_degredation_rate(poc, bgc, O₂, T)

    large_particle_degredation = λ * GOC
    degredation = λ * POC

    # grazing
    _, _, _, microzooplankton_grazing = specific_grazing(bgc.microzooplankton, P, D, Z, POC)
    _, _, _, mesozooplankton_grazing = specific_grazing(bgc.mesozooplankton, P, D, Z, POC)
    small_flux_feeding = specific_flux_feeding(bgc.mesozooplankton, POC, bgc.sinking_velocities.POC.w, grid)

    grazing = microzooplankton_grazing * Z + (mesozooplankton_grazing + small_flux_feeding) * M

    # aggregation
    _, Φ₁, _, Φ₃ = aggregation(bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, zₘₓₗ)
    dissolved_aggregation = Φ₁ + Φ₃
    
    aggregation_to_large = aggregation(poc, bgc, z, POC, GOC, zₘₓₗ)

    aggregation = dissolved_aggregation - aggregation_to_large

    return (grazing_waste 
            + nanophytoplankton_mortality + diatom_mortality + microzooplankton_mortality 
            + large_particle_degredation + aggregation 
            - grazing - degredation)
end

@inline function (poc::TwoCompartementParticulateOrganicMatter)(bgc, ::Val{:GOC}, 
                                                                x, y, z, t,
                                                                P, D, Z, M, 
                                                                PChl, DChl, PFe, DFe, DSi,
                                                                DOC, POC, GOC, 
                                                                SFe, BFe, PSi, 
                                                                NO₃, NH₄, PO₄, Fe, Si, 
                                                                CaCO₃, DIC, Alk, 
                                                                O₂, T, 
                                                                zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    grazing_waste = specific_non_assimilated_waste(bgc.mesozooplankton, bgc, P, D, Z, POC, GOC) * M

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, D, zₘₓₗ)

    nanophytoplankton_mortality = 0.5 * R_CaCO₃ * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality)

    diatom_linear_mortality, diatom_quadratic_mortality = mortality(bgc.diatoms, bgc, z, D, zₘₓₗ)

    diatom_mortality = 0.5 * diatom_linear_mortality + diatom_quadratic_mortality

    mesozooplankton_mortality = mortality(bgc.mesozooplankton, bgc, I, O₂, T)

    # degredation
    λ = specific_degredation_rate(poc, bgc, O₂, T)

    degredation = λ * GOC

    # grazing
    grazing = specific_flux_feeding(bgc.mesozooplankton, GOC, bgc.sinking_velocities.GOC.w, grid) * M

    # aggregation
    _, _, dissolved_aggregation = aggregation(bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, zₘₓₗ)
    
    small_particle_aggregation = aggregation(poc, bgc, z, POC, GOC, zₘₓₗ)

    aggregation = dissolved_aggregation + small_particle_aggregation

    # fecal pelet prodiction
    fecal_pelet_production = upper_trophic_fecal_product(bgc.mesozooplankton, M, T)

    return (grazing_waste 
            + nanophytoplankton_mortality + diatom_mortality + mesozooplankton_mortality 
            + aggregation + fecal_pelet_production
            - grazing 
            - degredation)
end
