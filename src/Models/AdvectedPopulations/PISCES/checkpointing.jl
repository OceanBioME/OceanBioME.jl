import Oceananigans: prognostic_state, restore_prognostic_state!

# the mixed layer depth and silicate climatology aren't computed by PISCES so may be set by the
# user (e.g. in a callback), while the other auxiliary fields are computed when the state is updated
prognostic_state(bgc::PISCES) =
    (mixed_layer_depth = field_state(bgc.mixed_layer_depth),
     silicate_climatology = field_state(bgc.silicate_climatology))

function restore_prognostic_state!(bgc::PISCES, from)
    restore_prognostic_state!(bgc.mixed_layer_depth, from.mixed_layer_depth)
    restore_prognostic_state!(bgc.silicate_climatology, from.silicate_climatology)
    return bgc
end

# `ConstantField`s, `FieldTimeSeries`, etc. have no state to checkpoint
field_state(field::Field) = prognostic_state(field)
field_state(field) = nothing
