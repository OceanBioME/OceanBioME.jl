using Adapt

import Adapt: adapt_structure

Adapt.adapt_structure(to, bnd::BiologyNutrientDetritus) =
    LOBSTER(adapt(to, bnd.nutrients),
            adapt(to, bnd.biology),
            adapt(to, bnd.detritus),
            adapt(to, bnd.carbonate_system),
            adapt(to, bnd.oxygen))

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
             