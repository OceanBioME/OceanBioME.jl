import Oceananigans: prognostic_state, restore_prognostic_state!

# the tracked fields are copied from the model when the state is updated
prognostic_state(sediment::BiogeochemicalSediment) =
    (fields = prognostic_state(sediment.fields),
     timestepper = prognostic_state(sediment.timestepper),
     clock = prognostic_state(sediment.clock))

function restore_prognostic_state!(sediment::BiogeochemicalSediment, from)
    restore_prognostic_state!(sediment.fields, from.fields)
    restore_prognostic_state!(sediment.timestepper, from.timestepper)
    restore_prognostic_state!(sediment.clock, from.clock)
    return sediment
end
