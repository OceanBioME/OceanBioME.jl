using Oceananigans.Architectures: CPU, on_architecture

import Oceananigans: prognostic_state, restore_prognostic_state!

# the tendencies are recomputed before the particles are stepped
prognostic_state(particles::BiogeochemicalParticles) =
    (x = on_architecture(CPU(), particles.x),
     y = on_architecture(CPU(), particles.y),
     z = on_architecture(CPU(), particles.z),
     fields = on_architecture(CPU(), particles.fields))

function restore_prognostic_state!(particles::BiogeochemicalParticles, from)
    restore_prognostic_state!(particles.x, from.x)
    restore_prognostic_state!(particles.y, from.y)
    restore_prognostic_state!(particles.z, from.z)
    restore_prognostic_state!(particles.fields, from.fields)
    return particles
end
