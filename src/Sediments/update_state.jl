using Oceananigans.BoundaryConditions: update_boundary_conditions!, fill_halo_regions!
using Oceananigans.TimeSteppers: time_step!

import Oceananigans.TimeSteppers: update_state!

function update_biogeochemical_state!(model, sediment_model::BiogeochemicalSediment)
    update_tracked_fields!(sediment_model, model)

    # the state can be updated without the model having been stepped (e.g. when a simulation is
    # initialised after picking up from a checkpoint, or by `set!`), so we only step the sediment
    # when the model clock has moved on since the sediment was last updated
    first_update = (model.clock.iteration == 0) && (model.clock.stage == 1)
    model_stepped = model.clock.time != sediment_model.clock.time

    if !first_update && model_stepped
        step_sediment!(sediment_model, model, model.timestepper)
    end

    sediment_model.clock.time = model.clock.time

    update_state!(sediment_model)

    return nothing
end

function update_state!(model::BiogeochemicalSediment, callbacks=[]; compute_tendencies = true)
    model_fields = fields(model)

    update_boundary_conditions!(model_fields, model)

    fill_halo_regions!(model_fields, model.clock, model_fields; async = true)

    compute_sediment_tendencies!(model)

    return nothing
end
