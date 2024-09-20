using Adapt

import Adapt: adapt_structure, adapt

# we can throw away all of the fields since they're delt with outside of the kernels
adapt_structure(to, bgc::PISCES) =
    PISCES(adapt(to, bgc.nanophytoplankton),
           adapt(to, bgc.diatoms),
           adapt(to, bgc.microzooplankton),
           adapt(to, bgc.mesozooplankton),
           adapt(to, bgc.dissolved_organic_matter),
           adapt(to, bgc.particulate_organic_matter),
           adapt(to, bgc.nitrogen),
           adapt(to, bgc.iron),
           adapt(to, bgc.silicate),
           adapt(to, bgc.oxygen),
           adapt(to, bgc.phosphate),
           adapt(to, bgc.calcite),
           adapt(to, bgc.carbon_system),
           adapt(to, bgc.first_anoxia_threshold),
           adapt(to, bgc.second_anoxia_threshold),
           adapt(to, bgc.nitrogen_redfield_ratio),
           adapt(to, bgc.phosphate_redfield_ratio),
           adapt(to, bgc.mixed_layer_shear),
           adapt(to, bgc.background),
           nothing, nothing, nothing, nothing, nothing,
           adapt(to, bgc.carbon_chemistry),
           nothing, nothing)

# adapting a bunch of numbers but maybe someone will want something else in there in the future
adapt_structure(to, phyto::Phytoplankton) = 
    Phytoplankton(adapt(to, phyto.growth_rate),
                  adapt(to, phyto.nutrient_limitation),
                  adapt(to, phyto.exudated_fracton).
                  adapt(to, phyto.blue_light_absorption),
                  adapt(to, phyto.green_light_absorption),
                  adapt(to, phyto.red_light_absorption),
                  adapt(to, phyto.mortality_half_saturation),
                  adapt(to, phyto.linear_mortality_rate),
                  adapt(to, phyto.base_quadratic_mortality),
                  adapt(to, phyto.maximum_quadratic_mortality),
                  adapt(to, phyto.minimum_chlorophyll_ratio),
                  adapt(to, phyto.maximum_chlorophyll_ratio),
                  adapt(to, phyto.maximum_iron_ratio),
                  adapt(to, phyto.silicate_half_saturation),
                  adapt(to, phyto.enhanced_silicate_half_saturation),
                  adapt(to, phyto.optimal_silicate_ratio))
