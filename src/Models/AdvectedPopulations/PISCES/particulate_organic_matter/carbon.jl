# these are just completly different to eachother so not going to try and define a generic function 

@inline function (bgc::TwoCompartementPOCPISCES)(i, j, k, grid, val_name::Val{:POC}, clock, fields, auxiliary_fields)
    # gains
    grazing_waste = small_non_assimilated_waste(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_mortality = small_mortality(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    zooplankton_mortality = small_mortality(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    _, Φ₁, _, Φ₃ = aggregation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    dissolved_aggregation = Φ₁ + Φ₃

    large_breakdown = degredation(bgc.particulate_organic_matter, Val(:GOC), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # losses
    grazing = total_grazing(bgc.zooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    aggregation_to_large = aggregation(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    small_breakdown = degredation(bgc.particulate_organic_matter, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (grazing_waste + phytoplankton_mortality + zooplankton_mortality + dissolved_aggregation + large_breakdown
            - grazing - aggregation_to_large - small_breakdown)
end

@inline function (bgc::TwoCompartementPOCPISCES)(i, j, k, grid, val_name::Val{:GOC}, clock, fields, auxiliary_fields)
    # gains
    grazing_waste = large_non_assimilated_waste(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_mortality = large_mortality(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    zooplankton_mortality = large_mortality(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    aggregation_to_large = aggregation(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic_feces = upper_trophic_fecal_production(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # losses
    grazing = total_grazing(bgc.zooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    large_breakdown = degredation(bgc.particulate_organic_matter, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (grazing_waste + phytoplankton_mortality + zooplankton_mortality + upper_trophic_feces
            - grazing  - large_breakdown)
end

@inline degredation(poc::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields, auxiliary_fields) = # for going to DOC
    degredation(poc::TwoCompartementCarbonIronParticles, Val(:POC), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline degredation(poc::TwoCompartementCarbonIronParticles, ::Val{:POC}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) = 
    @inbounds specific_degredation_rate(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields) * fields.POC[i, j, k]

@inline degredation(poc::TwoCompartementCarbonIronParticles, ::Val{:GOC}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) = 
    @inbounds specific_degredation_rate(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields) * fields.GOC[i, j, k]