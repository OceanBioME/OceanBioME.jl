using Adapt

import Adapt: adapt_structure

# we can throw away all of the fields since they're delt with outside of the kernels
Adapt.adapt_structure(to, bgc::PISCES) =
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
           adapt(to, bgc.background_shear),
           adapt(to, bgc.latitude),
           adapt(to, bgc.day_length),
           adapt(to, bgc.mixed_layer_depth),
           adapt(to, bgc.euphotic_depth), 
           adapt(to, bgc.silicate_climatology),
           adapt(to, bgc.mean_mixed_layer_vertical_diffusivity),
           adapt(to, bgc.mean_mixed_layer_light),
           adapt(to, bgc.carbon_chemistry),
           adapt(to, bgc.silicate_climatology), 
           adapt(to, bgc.sinking_velocities))
