using Adapt

import Adapt: adapt_structure

Adapt.adapt_structure(to, lobster::LOBSTER) =
    LOBSTER(adapt(to, lobster.nutrients),
            adapt(to, lobster.biology),
            adapt(to, lobster.detritus),
            adapt(to, lobster.carbonate_system),
            adapt(to, lobster.oxygen))

Adapt.adapt_structure(to, detritus::TwoParticleAndDissolved) =
    TwoParticleAndDissolved(adapt(to, detritus.remineralisation_inorganic_fraction),
                            adapt(to, detritus.small_reminerlisation_rate),
                            adapt(to, detritus.large_reminerlisation_rate),
                            adapt(to, detritus.dissolved_reminerlisation_rate),
                            adapt(to, detritus.small_solid_waste_fraction),
                            adapt(to, detritus.redfield_ratio),
                            adapt(to, detritus.small_particle_sinking_velocity),
                            adapt(to, detritus.large_particle_sinking_velocity))

Adapt.adapt_structure(to, detritus::VariableRedfieldDetritus) =
    VariableRedfieldDetritus(adapt(to, detritus.remineralisation_inorganic_fraction),
                             adapt(to, detritus.small_reminerlisation_rate),
                             adapt(to, detritus.large_reminerlisation_rate),
                             adapt(to, detritus.dissolved_reminerlisation_rate),
                             adapt(to, detritus.small_solid_waste_fraction),
                             adapt(to, detritus.small_particle_sinking_velocity),
                             adapt(to, detritus.large_particle_sinking_velocity))
