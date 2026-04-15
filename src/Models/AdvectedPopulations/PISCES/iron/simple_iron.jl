"""
    SimpleIron(; excess_scavenging_enhancement = 1000)

Parameterisation for iron evolution, not the "complex chemistry" model
of Aumount et al, 2015. Iron is scavenged (i.e. permanently removed from
the model) when the free iron concentration exceeds the ligand concentration
at a rate modified by `excess_scavenging_enhancement`.
"""
@kwdef struct SimpleIron{FT}
    excess_scavenging_enhancement :: FT = 1000.0 # unitless
     maximum_ligand_concentration :: FT = 0.6    # μmol Fe / m³
           dissolved_ligand_ratio :: FT = 0.09   # μmol Fe / mmol C
end

struct IronTendencyArgs{FT}
    Fe :: FT
    DOC :: FT
    T :: FT
    POC :: FT
    GOC :: FT
    SFe :: FT
    CaCO₃ :: FT
    PSi :: FT
    P :: FT
    PChl :: FT
    PFe :: FT
    D :: FT
    DChl :: FT
    DFe :: FT
    Z :: FT
    M :: FT
    NH₄ :: FT
    NO₃ :: FT
    PO₄ :: FT
    Si :: FT
    O₂ :: FT
    z :: FT
    zₘₓₗ :: FT
    zₑᵤ :: FT
    Si′ :: FT
    background_shear :: FT
    mixed_layer_shear :: FT
    sinking_flux :: FT
    sinking_iron_flux :: FT
    first_anoxia_threshold :: FT
    second_anoxia_threshold :: FT
    minimum_iron_scavenging_rate :: FT
    load_specific_iron_scavenging_rate :: FT
    base_breakdown_rate :: FT
    particle_temperature_sensitivity :: FT
    maximum_iron_ratio_in_bacteria :: FT
    iron_half_saturation_for_bacteria :: FT
    bacterial_iron_uptake_efficiency :: FT
    maximum_bacterial_growth_rate :: FT
    dissolved_organic_aggregation_parameter_1 :: FT
    dissolved_organic_aggregation_parameter_2 :: FT
    dissolved_organic_aggregation_parameter_3 :: FT
    dissolved_organic_aggregation_parameter_4 :: FT
    dissolved_organic_aggregation_parameter_5 :: FT
    microzooplankton_bacteria_concentration :: FT
    mesozooplankton_bacteria_concentration :: FT
    maximum_bacteria_concentration :: FT
    bacteria_concentration_depth_exponent :: FT
    doc_half_saturation_for_bacterial_activity :: FT
    nitrate_half_saturation_for_bacterial_activity :: FT
    ammonia_half_saturation_for_bacterial_activity :: FT
    phosphate_half_saturation_for_bacterial_activity :: FT
    iron_half_saturation_for_bacterial_activity :: FT
    nano_exudated_fraction :: FT
    nano_maximum_iron_ratio :: FT
    nano_half_saturation_for_iron_uptake :: FT
    nano_threshold_for_size_dependency :: FT
    nano_size_ratio :: FT
    nano_minimum_ammonium_half_saturation :: FT
    nano_minimum_nitrate_half_saturation :: FT
    nano_optimal_iron_quota :: FT
    nano_base_growth_rate :: FT
    nano_temperature_sensitivity :: FT
    diatom_exudated_fraction :: FT
    diatom_maximum_iron_ratio :: FT
    diatom_half_saturation_for_iron_uptake :: FT
    diatom_threshold_for_size_dependency :: FT
    diatom_size_ratio :: FT
    diatom_minimum_ammonium_half_saturation :: FT
    diatom_minimum_nitrate_half_saturation :: FT
    diatom_optimal_iron_quota :: FT
    diatom_base_growth_rate :: FT
    diatom_temperature_sensitivity :: FT
    microzooplankton_iron_ratio :: FT
    microzooplankton_non_assimilated_fraction :: FT
    microzooplankton_maximum_grazing_rate :: FT
    microzooplankton_temperature_sensitivity :: FT
    microzooplankton_preference_for_p :: FT
    microzooplankton_preference_for_d :: FT
    microzooplankton_preference_for_z :: FT
    microzooplankton_preference_for_poc :: FT
    microzooplankton_specific_food_threshold_concentration :: FT
    microzooplankton_grazing_half_saturation :: FT
    microzooplankton_food_threshold_concentration :: FT
    microzooplankton_minimum_growth_efficiency :: FT
    microzooplankton_maximum_flux_feeding_rate :: FT
    mesozooplankton_iron_ratio :: FT
    mesozooplankton_non_assimilated_fraction :: FT
    mesozooplankton_maximum_grazing_rate :: FT
    mesozooplankton_temperature_sensitivity :: FT
    mesozooplankton_preference_for_p :: FT
    mesozooplankton_preference_for_d :: FT
    mesozooplankton_preference_for_z :: FT
    mesozooplankton_preference_for_poc :: FT
    mesozooplankton_specific_food_threshold_concentration :: FT
    mesozooplankton_grazing_half_saturation :: FT
    mesozooplankton_food_threshold_concentration :: FT
    mesozooplankton_minimum_growth_efficiency :: FT
    mesozooplankton_maximum_flux_feeding_rate :: FT
    mesozooplankton_quadratic_mortality :: FT
end


@inline function iron_tendency(iron::SimpleIron, args::IronTendencyArgs)
    λ₀ = args.minimum_iron_scavenging_rate
    λ₁ = args.load_specific_iron_scavenging_rate
    λ₀ₚ = args.base_breakdown_rate
    bₚ = args.particle_temperature_sensitivity
    O₂ₘᵢₙ₁ = args.first_anoxia_threshold
    O₂ₘᵢₙ₂ = args.second_anoxia_threshold

    θᵦ = args.maximum_iron_ratio_in_bacteria
    Kᵦ = args.iron_half_saturation_for_bacteria
    κᵦ = args.bacterial_iron_uptake_efficiency
    μ₀ᵦ = args.maximum_bacterial_growth_rate

    a₁ = args.dissolved_organic_aggregation_parameter_1
    a₂ = args.dissolved_organic_aggregation_parameter_2
    a₃ = args.dissolved_organic_aggregation_parameter_3
    a₄ = args.dissolved_organic_aggregation_parameter_4
    a₅ = args.dissolved_organic_aggregation_parameter_5

    bZ = args.microzooplankton_bacteria_concentration
    bM = args.mesozooplankton_bacteria_concentration
    Bₘₐₓ = args.maximum_bacteria_concentration
    a = args.bacteria_concentration_depth_exponent
    K_DOC = args.doc_half_saturation_for_bacterial_activity
    K_NO₃ = args.nitrate_half_saturation_for_bacterial_activity
    K_NH₄ = args.ammonia_half_saturation_for_bacterial_activity
    K_PO₄ = args.phosphate_half_saturation_for_bacterial_activity
    K_Fe = args.iron_half_saturation_for_bacterial_activity

    λFe = iron_scavenging_rate(λ₀,
                               λ₁,
                               args.POC,
                               args.GOC,
                               args.CaCO₃,
                               args.PSi)
    Fe′ = free_iron(iron, args.Fe, args.DOC, args.T)

    food_availability = (; P = args.P,
                           D = args.D,
                           POC = args.POC,
                           Z = args.Z)
    iron_availability = (; P = args.PFe / (args.P + eps(zero(args.P))),
                           D = args.DFe / (args.D + eps(zero(args.D))),
                           POC = args.SFe / (args.POC + eps(zero(args.POC))),
                           Z = args.microzooplankton_iron_ratio)

    Bact = bacteria_concentration(bZ,
                                  bM,
                                  Bₘₐₓ,
                                  a,
                                  args.z,
                                  args.zₘₓₗ,
                                  args.zₑᵤ,
                                  args.Z,
                                  args.M)
    LBact = bacteria_activity(K_DOC,
                              K_NO₃,
                              K_NH₄,
                              K_PO₄,
                              K_Fe,
                              args.NH₄,
                              args.NO₃,
                              args.PO₄,
                              args.Fe,
                              args.DOC)

    small_particle_iron_remineralisation = degradation(Val(:SFe),
                                                       specific_degradation_rate(λ₀ₚ,
                                                                                 bₚ,
                                                                                 args.O₂,
                                                                                 args.T,
                                                                                 O₂ₘᵢₙ₁,
                                                                                 O₂ₘᵢₙ₂),
                                                       args.SFe)
    grazing_waste = non_assimilated_iron(args.microzooplankton_iron_ratio,
                                         args.microzooplankton_non_assimilated_fraction,
                                         args.microzooplankton_maximum_grazing_rate,
                                         args.microzooplankton_temperature_sensitivity,
                                         args.microzooplankton_preference_for_p,
                                         args.microzooplankton_preference_for_d,
                                         args.microzooplankton_preference_for_z,
                                         args.microzooplankton_preference_for_poc,
                                         args.microzooplankton_specific_food_threshold_concentration,
                                         args.microzooplankton_grazing_half_saturation,
                                         args.microzooplankton_food_threshold_concentration,
                                         args.microzooplankton_minimum_growth_efficiency,
                                         args.microzooplankton_maximum_flux_feeding_rate,
                                         args.mesozooplankton_iron_ratio,
                                         args.mesozooplankton_non_assimilated_fraction,
                                         args.mesozooplankton_maximum_grazing_rate,
                                         args.mesozooplankton_temperature_sensitivity,
                                         args.mesozooplankton_preference_for_p,
                                         args.mesozooplankton_preference_for_d,
                                         args.mesozooplankton_preference_for_z,
                                         args.mesozooplankton_preference_for_poc,
                                         args.mesozooplankton_specific_food_threshold_concentration,
                                         args.mesozooplankton_grazing_half_saturation,
                                         args.mesozooplankton_food_threshold_concentration,
                                         args.mesozooplankton_minimum_growth_efficiency,
                                         args.mesozooplankton_maximum_flux_feeding_rate,
                                         args.T,
                                         args.Z,
                                         args.M,
                                         food_availability,
                                         iron_availability,
                                         args.sinking_flux,
                                         args.sinking_iron_flux)
    upper_trophic_waste = upper_trophic_dissolved_iron(args.mesozooplankton_minimum_growth_efficiency,
                                                       args.mesozooplankton_non_assimilated_fraction,
                                                       args.mesozooplankton_iron_ratio,
                                                       args.mesozooplankton_temperature_sensitivity,
                                                       args.mesozooplankton_quadratic_mortality,
                                                       args.T,
                                                       args.M)
    phytoplankton_iron_uptake = uptake(args.nano_exudated_fraction,
                                       args.nano_maximum_iron_ratio,
                                       args.nano_half_saturation_for_iron_uptake,
                                       args.nano_threshold_for_size_dependency,
                                       args.nano_size_ratio,
                                       args.nano_minimum_ammonium_half_saturation,
                                       args.nano_minimum_nitrate_half_saturation,
                                       args.nano_optimal_iron_quota,
                                       args.nano_base_growth_rate,
                                       args.nano_temperature_sensitivity,
                                       args.diatom_exudated_fraction,
                                       args.diatom_maximum_iron_ratio,
                                       args.diatom_half_saturation_for_iron_uptake,
                                       args.diatom_threshold_for_size_dependency,
                                       args.diatom_size_ratio,
                                       args.diatom_minimum_ammonium_half_saturation,
                                       args.diatom_minimum_nitrate_half_saturation,
                                       args.diatom_optimal_iron_quota,
                                       args.diatom_base_growth_rate,
                                       args.diatom_temperature_sensitivity,
                                       Val(:Fe),
                                       args.T,
                                       args.Fe,
                                       args.NO₃,
                                       args.NH₄,
                                       args.PO₄,
                                       args.Si,
                                       args.Si′,
                                       args.P,
                                       args.PChl,
                                       args.PFe,
                                       args.D,
                                       args.DChl,
                                       args.DFe)
    colloidal_aggregation, = aggregation_of_colloidal_iron(a₁,
                                                           a₂,
                                                           a₃,
                                                           a₄,
                                                           a₅,
                                                           args.background_shear,
                                                           args.mixed_layer_shear,
                                                           args.z,
                                                           args.zₘₓₗ,
                                                           args.Fe,
                                                           Fe′,
                                                           args.DOC,
                                                           args.POC,
                                                           args.GOC)
    scavenging = iron_scavenging(λFe, args.POC + args.GOC, Fe′)
    bacterial_uptake = bacterial_iron_uptake(μ₀ᵦ,
                                             bₚ,
                                             θᵦ,
                                             Kᵦ,
                                             κᵦ,
                                             args.T,
                                             args.Fe,
                                             Bact,
                                             LBact)

    ligand_aggregation_loss = ligand_aggregation(iron, args.Fe, args.DOC, args.T, λFe)

    return (small_particle_iron_remineralisation + grazing_waste + upper_trophic_waste -
            phytoplankton_iron_uptake - ligand_aggregation_loss - colloidal_aggregation -
            scavenging - bacterial_uptake)
end

required_biogeochemical_tracers(::SimpleIron) = tuple(:Fe)

const SimpleIronPISCES = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:SimpleIron}

@inline function iron_tendency(bgc::SimpleIronPISCES, i, j, k, grid, clock, fields, auxiliary_fields)
    Fe = @inbounds fields.Fe[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]
    T = @inbounds fields.T[i, j, k]
    O₂ = @inbounds fields.O₂[i, j, k]

    POC = @inbounds fields.POC[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]
    SFe = @inbounds fields.SFe[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]
    PSi = @inbounds fields.PSi[i, j, k]

    P = @inbounds fields.P[i, j, k]
    PChl = @inbounds fields.PChl[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D = @inbounds fields.D[i, j, k]
    DChl = @inbounds fields.DChl[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    Z = @inbounds fields.Z[i, j, k]
    M = @inbounds fields.M[i, j, k]
    NH₄ = @inbounds fields.NH₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    PO₄ = @inbounds fields.PO₄[i, j, k]
    Si = @inbounds fields.Si[i, j, k]

    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]
    zₑᵤ = @inbounds auxiliary_fields.zₑᵤ[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())
    Si′ = @inbounds bgc.silicate_climatology[i, j, k]

    sinking_flux = edible_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)
    sinking_iron_flux = edible_iron_flux_rate(bgc.particulate_organic_matter, i, j, k, grid, fields, auxiliary_fields)

    pom = bgc.particulate_organic_matter
    (a₁,
     a₂,
     a₃,
     a₄,
     a₅) = bgc.dissolved_organic_matter.aggregation_parameters

    nano = bgc.phytoplankton.nano
    diatoms = bgc.phytoplankton.diatoms
    zoo = bgc.zooplankton
    micro = zoo.micro
    meso = zoo.meso

    args = IronTendencyArgs(
        Fe,
        DOC,
        T,
        POC,
        GOC,
        SFe,
        CaCO₃,
        PSi,
        P,
        PChl,
        PFe,
        D,
        DChl,
        DFe,
        Z,
        M,
        NH₄,
        NO₃,
        PO₄,
        Si,
        O₂,
        z,
        zₘₓₗ,
        zₑᵤ,
        Si′,
        bgc.background_shear,
        bgc.mixed_layer_shear,
        sinking_flux,
        sinking_iron_flux,
        bgc.first_anoxia_threshold,
        bgc.second_anoxia_threshold,
        pom.minimum_iron_scavenging_rate,
        pom.load_specific_iron_scavenging_rate,
        pom.base_breakdown_rate,
        pom.temperature_sensitivity,
        pom.maximum_iron_ratio_in_bacteria,
        pom.iron_half_saturation_for_bacteria,
        pom.bacterial_iron_uptake_efficiency,
        pom.maximum_bacterial_growth_rate,
        a₁,
        a₂,
        a₃,
        a₄,
        a₅,
        zoo.microzooplankton_bacteria_concentration,
        zoo.mesozooplankton_bacteria_concentration,
        zoo.maximum_bacteria_concentration,
        zoo.bacteria_concentration_depth_exponent,
        zoo.doc_half_saturation_for_bacterial_activity,
        zoo.nitrate_half_saturation_for_bacterial_activity,
        zoo.ammonia_half_saturation_for_bacterial_activity,
        zoo.phosphate_half_saturation_for_bacterial_activity,
        zoo.iron_half_saturation_for_bacterial_activity,
        nano.exudated_fraction,
        nano.maximum_iron_ratio,
        nano.half_saturation_for_iron_uptake,
        nano.threshold_for_size_dependency,
        nano.size_ratio,
        nano.nutrient_limitation.minimum_ammonium_half_saturation,
        nano.nutrient_limitation.minimum_nitrate_half_saturation,
        nano.nutrient_limitation.optimal_iron_quota,
        nano.growth_rate.base_growth_rate,
        nano.growth_rate.temperature_sensitivity,
        diatoms.exudated_fraction,
        diatoms.maximum_iron_ratio,
        diatoms.half_saturation_for_iron_uptake,
        diatoms.threshold_for_size_dependency,
        diatoms.size_ratio,
        diatoms.nutrient_limitation.minimum_ammonium_half_saturation,
        diatoms.nutrient_limitation.minimum_nitrate_half_saturation,
        diatoms.nutrient_limitation.optimal_iron_quota,
        diatoms.growth_rate.base_growth_rate,
        diatoms.growth_rate.temperature_sensitivity,
        micro.iron_ratio,
        micro.non_assimilated_fraction,
        micro.maximum_grazing_rate,
        micro.temperature_sensitivity,
        micro.food_preferences.P,
        micro.food_preferences.D,
        micro.food_preferences.Z,
        micro.food_preferences.POC,
        micro.specific_food_threshold_concentration,
        micro.grazing_half_saturation,
        micro.food_threshold_concentration,
        micro.minimum_growth_efficiency,
        micro.maximum_flux_feeding_rate,
        meso.iron_ratio,
        meso.non_assimilated_fraction,
        meso.maximum_grazing_rate,
        meso.temperature_sensitivity,
        meso.food_preferences.P,
        meso.food_preferences.D,
        meso.food_preferences.Z,
        meso.food_preferences.POC,
        meso.specific_food_threshold_concentration,
        meso.grazing_half_saturation,
        meso.food_threshold_concentration,
        meso.minimum_growth_efficiency,
        meso.maximum_flux_feeding_rate,
        meso.quadratic_mortality,
    )

    return iron_tendency(bgc.iron, args)
end

@inline function (bgc::SimpleIronPISCES)(i, j, k, grid, ::Val{:Fe}, clock, fields, auxiliary_fields)
    return iron_tendency(bgc, i, j, k, grid, clock, fields, auxiliary_fields)
end
