using Adapt

import Adapt: adapt_structure

Adapt.adapt_structure(to, bgc::PISCES) =
    PISCES(adapt(to, bgc.phytoplankton),
           adapt(to, bgc.zooplankton),
           adapt(to, bgc.dissolved_organic_matter),
           adapt(to, bgc.particulate_organic_matter),
           adapt(to, bgc.nitrogen),
           adapt(to, bgc.iron),
           adapt(to, bgc.silicate),
           adapt(to, bgc.oxygen),
           adapt(to, bgc.phosphate),
           adapt(to, bgc.inorganic_carbon),
           adapt(to, bgc.first_anoxia_threshold),
           adapt(to, bgc.second_anoxia_threshold),
           adapt(to, bgc.nitrogen_redfield_ratio),
           adapt(to, bgc.phosphate_redfield_ratio),
           adapt(to, bgc.mixed_layer_shear),
           adapt(to, bgc.background_shear),
           adapt(to, bgc.latitude),
           adapt(to, bgc.day_length),
           adapt(to, bgc.mixed_layer_depth),
           adapt(to, bgc.euphotic_depth),
           adapt(to, bgc.silicate_climatology),
           adapt(to, bgc.mean_mixed_layer_vertical_diffusivity),
           adapt(to, bgc.mean_mixed_layer_light),
           adapt(to, bgc.carbon_chemistry),
           adapt(to, bgc.calcite_saturation),
           adapt(to, bgc.sinking_velocities))

Adapt.adapt_structure(to, zoo::MicroAndMeso) =
    MicroAndMeso(adapt(to, zoo.micro),
                 adapt(to, zoo.meso),
                 adapt(to, zoo.microzooplankton_bacteria_concentration),
                 adapt(to, zoo.mesozooplankton_bacteria_concentration),
                 adapt(to, zoo.maximum_bacteria_concentration),
                 adapt(to, zoo.bacteria_concentration_depth_exponent),
                 adapt(to, zoo.doc_half_saturation_for_bacterial_activity),
                 adapt(to, zoo.nitrate_half_saturation_for_bacterial_activity),
                 adapt(to, zoo.ammonia_half_saturation_for_bacterial_activity),
                 adapt(to, zoo.phosphate_half_saturation_for_bacterial_activity),
                 adapt(to, zoo.iron_half_saturation_for_bacterial_activity))

Adapt.adapt_structure(to, zoo::QualityDependantZooplankton) =
    QualityDependantZooplankton(adapt(to, zoo.temperature_sensetivity),
                                adapt(to, zoo.maximum_grazing_rate),
                                adapt(to, zoo.food_preferences), # the only one that isn't already bits
                                adapt(to, zoo.food_threshold_concentration),
                                adapt(to, zoo.specific_food_thresehold_concentration),
                                adapt(to, zoo.grazing_half_saturation),
                                adapt(to, zoo.maximum_flux_feeding_rate),
                                adapt(to, zoo.iron_ratio),
                                adapt(to, zoo.minimum_growth_efficiency),
                                adapt(to, zoo.non_assililated_fraction),
                                adapt(to, zoo.mortality_half_saturation),
                                adapt(to, zoo.quadratic_mortality),
                                adapt(to, zoo.linear_mortality),
                                adapt(to, zoo.dissolved_excretion_fraction),
                                adapt(to, zoo.undissolved_calcite_fraction))
