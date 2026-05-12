using Adapt

import Adapt: adapt_structure

Adapt.adapt_structure(to, npd::NutrientsPlanktonDetritus) =
    NutrientsPlanktonDetritus(adapt(to, npd.nutrients),
                              adapt(to, npd.plankton),
                              adapt(to, npd.detritus),
                              adapt(to, npd.carbonate_system),
                              adapt(to, npd.oxygen))

Adapt.adapt_structure(to, detritus::TwoParticleAndDissolved) =
    TwoParticleAndDissolved(adapt(to, detritus.remineralisation_inorganic_fraction),
                            adapt(to, detritus.small_remineralisation_rate),
                            adapt(to, detritus.large_remineralisation_rate),
                            adapt(to, detritus.dissolved_remineralisation_rate),
                            adapt(to, detritus.small_solid_waste_fraction),
                            adapt(to, detritus.redfield_ratio),
                            adapt(to, detritus.small_particle_sinking_velocity),
                            adapt(to, detritus.large_particle_sinking_velocity))

Adapt.adapt_structure(to, detritus::VariableRedfieldDetritus) =
    VariableRedfieldDetritus(adapt(to, detritus.remineralisation_inorganic_fraction),
                             adapt(to, detritus.small_remineralisation_rate),
                             adapt(to, detritus.large_remineralisation_rate),
                             adapt(to, detritus.dissolved_remineralisation_rate),
                             adapt(to, detritus.small_solid_waste_fraction),
                             adapt(to, detritus.small_particle_sinking_velocity),
                             adapt(to, detritus.large_particle_sinking_velocity))

Adapt.adapt_structure(to, detritus::Detritus) =
    Detritus(adapt(to, detritus.remineralisation_rate),
             adapt(to, detritus.small_particle_fraction),
             adapt(to, detritus.redfield_ratio),
             adapt(to, detritus.sinking_speeds))
             
Adapt.adapt_structure(to, pz::PhytoZoo) =
    PhytoZoo(adapt(to, pz.nitrate_half_saturation),
             adapt(to, pz.ammonia_half_saturation),
             adapt(to, pz.iron_half_saturation),
             adapt(to, pz.nitrate_ammonia_inhibition),
             adapt(to, pz.light_half_saturation),
             adapt(to, pz.phytoplankton_maximum_growth_rate),
             adapt(to, pz.iron_ratio),
             adapt(to, pz.phytoplankton_exudation_fraction),
             adapt(to, pz.ammonia_fraction_of_exudate),
             adapt(to, pz.light_limitation),
             adapt(to, pz.temperature_coefficient),
             adapt(to, pz.phytoplankton_mortality_rate),
             adapt(to, pz.zooplankton_mortality_rate),
             adapt(to, pz.zooplankton_excretion_rate),
             adapt(to, pz.phytoplankton_mortality_formulation),
             adapt(to, pz.phytoplankton_solid_waste_fraction),
             adapt(to, pz.excretion_inorganic_fraction),
             adapt(to, pz.preference_for_phytoplankton),
             adapt(to, pz.maximum_grazing_rate),
             adapt(to, pz.grazing_half_saturation),
             adapt(to, pz.zooplankton_assimilation_fraction),
             adapt(to, pz.grazing_concentration_formulation),
             adapt(to, pz.zooplankton_calcite_dissolution),
             adapt(to, pz.redfield_ratio),
             adapt(to, pz.carbon_calcite_ratio),
             adapt(to, pz.zooplankton_gut_calcite_dissolution),
             adapt(to, pz.phytoplankton_chlorophyll_ratio),
             adapt(to, pz.phytoplankton_sinking_velocity),
             adapt(to, pz.zooplankton_sinking_velocity))