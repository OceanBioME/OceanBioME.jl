
@inline iron_ratio(iron_inventory, carbon_inventory) = iron_inventory / (carbon_inventory + eps(0.0))

@inline function small_particulate_iron_tendency(POC,
                                                 SFe,
                                                 grazing_waste,
                                                 phytoplankton_mortality,
                                                 zooplankton_mortality,
                                                 large_breakdown,
                                                 scavenging,
                                                 κ,
                                                 BactFe,
                                                 colloidal_aggregation,
                                                 grazing,
                                                 aggregation_to_large,
                                                 small_breakdown)

    θ = iron_ratio(SFe, POC)

    bacterial_assimilation = κ * BactFe
    grazing_iron = grazing * θ
    aggregation_to_large_iron = aggregation_to_large * θ

    return (grazing_waste + phytoplankton_mortality + zooplankton_mortality
            + large_breakdown + scavenging + bacterial_assimilation + colloidal_aggregation
            - grazing_iron - aggregation_to_large_iron - small_breakdown)
end

@inline function large_particulate_iron_tendency(POC,
                                                 SFe,
                                                 GOC,
                                                 BFe,
                                                 grazing_waste,
                                                 phytoplankton_mortality,
                                                 zooplankton_mortality,
                                                 upper_trophic_feces,
                                                 scavenging,
                                                 κ,
                                                 BactFe,
                                                 colloidal_aggregation,
                                                 aggregation_to_large,
                                                 grazing,
                                                 large_breakdown)
    
    θS = iron_ratio(SFe, POC)
    θB = iron_ratio(BFe, GOC)

    bacterial_assimilation = κ * BactFe
    grazing_iron = grazing * θB
    aggregation_to_large_iron = aggregation_to_large * θS

    return (grazing_waste + phytoplankton_mortality + zooplankton_mortality + upper_trophic_feces
            + scavenging + bacterial_assimilation + colloidal_aggregation + aggregation_to_large_iron
            - grazing_iron - large_breakdown)
end

@inline function (bgc::TwoCompartmentPOCPISCES)(i, j, k, grid, val_name::Val{:SFe}, clock, fields, auxiliary_fields)
    POC = @inbounds fields.POC[i, j, k]
    SFe = @inbounds fields.SFe[i, j, k]

    # gains
    grazing_waste = 
        small_non_assimilated_iron_waste(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_mortality = 
        small_mortality_iron(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    zooplankton_mortality = 
        small_mortality_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    large_breakdown = 
        degradation(bgc.particulate_organic_matter, Val(:BFe), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    λFe = iron_scavenging_rate(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    Fe′ = free_iron(bgc.iron, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    scavenging = iron_scavenging(λFe, POC, Fe′)

    κ = bgc.particulate_organic_matter.small_fraction_of_bacterially_consumed_iron

    BactFe = bacterial_iron_uptake(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    _, colloidal_aggregation = aggregation_of_colloidal_iron(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # losses
    grazing = total_grazing(bgc.zooplankton, Val(:POC), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    aggregation_to_large = aggregation(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    small_breakdown = degradation(bgc.particulate_organic_matter, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return small_particulate_iron_tendency(POC,
                                           SFe,
                                           grazing_waste,
                                           phytoplankton_mortality,
                                           zooplankton_mortality,
                                           large_breakdown,
                                           scavenging,
                                           κ,
                                           BactFe,
                                           colloidal_aggregation,
                                           grazing,
                                           aggregation_to_large,
                                           small_breakdown)
end

@inline function (bgc::TwoCompartmentPOCPISCES)(i, j, k, grid, val_name::Val{:BFe}, clock, fields, auxiliary_fields)
    POC = @inbounds fields.POC[i, j, k]
    SFe = @inbounds fields.SFe[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]
    BFe = @inbounds fields.BFe[i, j, k]

    # gains
    grazing_waste = large_non_assimilated_iron_waste(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_mortality = large_mortality_iron(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    zooplankton_mortality = large_mortality_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    aggregation_to_large = aggregation(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic_feces = upper_trophic_fecal_iron_production(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    λFe = iron_scavenging_rate(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    Fe′ = free_iron(bgc.iron, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    scavenging = iron_scavenging(λFe, GOC, Fe′)

    κ = bgc.particulate_organic_matter.large_fraction_of_bacterially_consumed_iron

    BactFe = bacterial_iron_uptake(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    _, _, colloidal_aggregation = aggregation_of_colloidal_iron(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # losses
    grazing = total_grazing(bgc.zooplankton, Val(:GOC), i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    large_breakdown = degradation(bgc.particulate_organic_matter, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return large_particulate_iron_tendency(POC,
                                           SFe,
                                           GOC,
                                           BFe,
                                           grazing_waste,
                                           phytoplankton_mortality,
                                           zooplankton_mortality,
                                           upper_trophic_feces,
                                           scavenging,
                                           κ,
                                           BactFe,
                                           colloidal_aggregation,
                                           aggregation_to_large,
                                           grazing,
                                           large_breakdown)
end

@inline degradation(::Val{:SFe}, degradation_rate, SFe) =
    degradation_rate * SFe

@inline degradation(::Val{:BFe}, degradation_rate, BFe) =
    degradation_rate * BFe

@inline degradation(poc::TwoCompartmentCarbonIronParticles, ::Val{:SFe}, degradation_rate, SFe) =
    degradation(Val(:SFe), degradation_rate, SFe)

@inline degradation(poc::TwoCompartmentCarbonIronParticles, ::Val{:BFe}, degradation_rate, BFe) =
    degradation(Val(:BFe), degradation_rate, BFe)

@inline degradation(poc::TwoCompartmentCarbonIronParticles, ::Val{:SFe}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    @inbounds degradation(poc, Val(:SFe), specific_degradation_rate(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields), fields.SFe[i, j, k])

@inline degradation(poc::TwoCompartmentCarbonIronParticles, ::Val{:BFe}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    @inbounds degradation(poc, Val(:BFe), specific_degradation_rate(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields), fields.BFe[i, j, k])

@inline iron_scavenging_rate(λ₀,
                                     λ₁,
                                     POC,
                                     GOC,
                                     CaCO₃,
                                     PSi) =
    λ₀ + λ₁ * (POC + GOC + CaCO₃ + PSi)

@inline function iron_scavenging_rate(pom::TwoCompartmentCarbonIronParticles, POC, GOC, CaCO₃, PSi)
    return iron_scavenging_rate(pom.minimum_iron_scavenging_rate,
                                pom.load_specific_iron_scavenging_rate,
                                POC,
                                GOC,
                                CaCO₃,
                                PSi)
end

@inline function iron_scavenging_rate(pom::TwoCompartmentCarbonIronParticles, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    POC = @inbounds fields.POC[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]
    PSi = @inbounds fields.PSi[i, j, k]

    return iron_scavenging_rate(pom, POC, GOC, CaCO₃, PSi)
end

@inline function bacterial_iron_uptake(μ₀,
                                       b,
                                       θ,
                                       K,
                                       κ,
                                       T,
                                       Fe,
                                       Bact,
                                       LBact)
    μ = μ₀ * b^T

    return μ * LBact * θ * Fe / (Fe + K) * Bact * κ
end

@inline function bacterial_iron_uptake(poc::TwoCompartmentCarbonIronParticles, T, Fe, Bact, LBact)
    return bacterial_iron_uptake(poc.maximum_bacterial_growth_rate,
                                 poc.temperature_sensitivity,
                                 poc.maximum_iron_ratio_in_bacteria,
                                 poc.iron_half_saturation_for_bacteria,
                                 poc.bacterial_iron_uptake_efficiency,
                                 T,
                                 Fe,
                                 Bact,
                                 LBact)
end

@inline function bacterial_iron_uptake(poc::TwoCompartmentCarbonIronParticles, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    T = @inbounds fields.T[i, j, k]
    Fe = @inbounds fields.Fe[i, j, k]

    Bact = bacteria_concentration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    LBact = bacteria_activity(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return bacterial_iron_uptake(poc, T, Fe, Bact, LBact)
end

@inline iron_scavenging(λFe, particle_load, Fe′) = λFe * particle_load * Fe′

@inline function iron_scavenging(poc::TwoCompartmentCarbonIronParticles, POC, GOC, CaCO₃, PSi, Fe′)
    λFe = iron_scavenging_rate(poc, POC, GOC, CaCO₃, PSi)

    return iron_scavenging(λFe, POC + GOC, Fe′)
end

