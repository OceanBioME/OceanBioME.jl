using Oceananigans.BoundaryConditions: update_boundary_conditions!, fill_halo_regions!
using Oceananigans.TimeSteppers: time_step!

import Oceananigans.TimeSteppers: update_state!

function update_biogeochemical_state!(model, sediment_model::BiogeochemicalSediment)
    update_tracked_fields!(sediment_model, model)

    priming_call = (model.clock.iteration == 0) & (model.clock.stage == 1)

    if !priming_call
        step_sediment!(sediment_model, model, model.timestepper)
    end

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
