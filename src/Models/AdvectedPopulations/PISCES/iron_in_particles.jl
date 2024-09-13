@inline function (doc::TwoCompartementParticulateOrganicMatter)(bgc, ::Val{:SFe}, 
                                                                x, y, z, t,
                                                                P, D, Z, M, 
                                                                PChl, DChl, PFe, DFe, DSi,
                                                                DOC, POC, GOC, 
                                                                SFe, BFe, PSi, 
                                                                NO₃, NH₄, PO₄, Fe, Si, 
                                                                CaCO₃, DIC, Alk, 
                                                                O₂, T, 
                                                                zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    grazing_waste = specific_non_assimilated_iron_waste(bgc.microzooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, BFe) * Z

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, D, zₘₓₗ)

    nanophytoplankton_mortality = (1 - 0.5 * R_CaCO₃) * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality) * PFe / P

    diatom_linear_mortality, = mortality(bgc.diatoms, bgc, z, D, zₘₓₗ)

    diatom_mortality = 0.5 * diatom_linear_mortality * DFe / D

    microzooplankton_mortality = mortality(bgc.microzooplankton, bgc, I, O₂, T) * bgc.microzooplankton.iron_ratio

    # degredation
    λ = specific_degredation_rate(doc, bgc, O₂, T)

    large_particle_degredation = λ * BFe
    degredation = λ * SFe

    # grazing
    _, _, _, microzooplankton_grazing = specific_grazing(bgc.microzooplankton, P, D, Z, POC)
    _, _, _, mesozooplankton_grazing = specific_grazing(bgc.mesozooplankton, P, D, Z, POC)
    small_flux_feeding = specific_flux_feeding(bgc.mesozooplankton, POC, bgc.sinking_velocities.POC.w, grid)

    grazing = (microzooplankton_grazing * Z + (mesozooplankton_grazing + small_flux_feeding) * M) * SFe / POC

    # aggregation
    
    aggregation_to_large = aggregation(doc, bgc, z, POC, GOC, zₘₓₗ)

    aggregation = aggregation_to_large * SFe / POC

    # scavenging
    λ₀ = doc.minimum_iron_scavenging_rate
    λ₁ = doc.load_specific_iron_scavenging_rate
    λ₂ = doc.dust_specific_iron_scavenging_rate

    λFe = iron_scavenging_rate(doc, POC, GOC, CaCO₃, PSi, dust)
    
    Fe′ = free_iron(bgc.iron, Fe, DOC, T)

    scavenging = λFe * POC * Fe′

    # bacterial uptake of dissolved iron
    κ = small_fraction_of_bacterially_consumed_iron

    BactFe = bacterial_iron_uptake(bgc.dissolved_organic_matter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    bacterial_assimilation = κ * BactFe

    # colloidal iron aggregation
    _, colloidal_aggregation = aggregation_of_colloidal_iron(bgc.iron, bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, Fe, T, zₘₓₗ)

    return (grazing_waste 
            + nanophytoplankton_mortality + diatom_mortality + microzooplankton_mortality 
            + large_particle_degredation + scavenging + bacterial_assimilation + colloidal_aggregation
            - aggregation 
            - grazing - degredation)
end


@inline function (poc::TwoCompartementParticulateOrganicMatter)(bgc, ::Val{:BFe}, 
                                                                x, y, z, t,
                                                                P, D, Z, M, 
                                                                PChl, DChl, PFe, DFe, DSi,
                                                                DOC, POC, GOC, 
                                                                SFe, BFe, PSi, 
                                                                NO₃, NH₄, PO₄, Fe, Si, 
                                                                CaCO₃, DIC, Alk, 
                                                                O₂, T, 
                                                                zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    grazing_waste = specific_non_assimilated_iron_waste(bgc.mesozooplankton, P, D, PFe, DFe, Z, POC, GOC, SFe, BFe) * M

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, D, zₘₓₗ)

    nanophytoplankton_mortality = 0.5 * R_CaCO₃ * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality) * PFe / P

    diatom_linear_mortality, diatom_quadratic_mortality = mortality(bgc.diatoms, bgc, z, D, zₘₓₗ)

    diatom_mortality = (0.5 * diatom_linear_mortality + diatom_quadratic_mortality) * DFe / D

    mesozooplankton_mortality = mortality(bgc.mesozooplankton, bgc, I, O₂, T) * bgc.mesozooplankton.iron_ratio

    # degredation
    λ = specific_degredation_rate(poc, bgc, O₂, T)

    degredation = λ * BFe

    # grazing
    grazing = specific_flux_feeding(bgc.mesozooplankton, GOC, bgc.sinking_velocities.GOC.w, grid) * M * SFe / POC

    # aggregation    
    small_particle_aggregation = aggregation(poc, bgc, z, POC, GOC, zₘₓₗ) 

    aggregation = small_particle_aggregation * SFe / POC

    # fecal pelet prodiction
    fecal_pelet_production = upper_trophic_fecal_product(bgc.mesozooplankton, M, T) * bgc.mesozooplankton.iron_ratio

    # scavenging
    λ₀ = doc.minimum_iron_scavenging_rate
    λ₁ = doc.load_specific_iron_scavenging_rate
    λ₂ = doc.dust_specific_iron_scavenging_rate

    λFe = λ₀ + λ₁ * (POC + GOC + CaCO₃ + PSi) + λ₂ * dust
    
    Fe′ = free_iron(bgc.iron, Fe, DOC, T)

    scavenging = λFe * GOC * Fe′

    # bacterial uptake of dissolved iron
    κ = large_fraction_of_bacterially_consumed_iron

    BactFe = bacterial_iron_uptake(bgc.dissolved_organic_matter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    bacterial_assimilation = κ * BactFe

    # colloidal iron aggregation
    _, _, colloidal_aggregation = aggregation_of_colloidal_iron(bgc.iron, bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, Fe, T, zₘₓₗ)

    return (grazing_waste
            + nanophytoplankton_mortality + diatom_mortality + mesozooplankton_mortality 
            + aggregation + fecal_pelet_production + scavenging + bacterial_assimilation + colloidal_aggregation
            - grazing - degredation)
end

@inline function iron_scavenging_rate(doc, POC, GOC, CaCO₃, PSi, dust)
    λ₀ = doc.minimum_iron_scavenging_rate
    λ₁ = doc.load_specific_iron_scavenging_rate
    λ₂ = doc.dust_specific_iron_scavenging_rate

    return λ₀ + λ₁ * (POC + GOC + CaCO₃ + PSi) + λ₂ * dust
end