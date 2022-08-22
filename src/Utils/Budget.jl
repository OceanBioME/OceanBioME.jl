module Budget
using OceanBioME
using Oceananigans.Operators: Vᶜᶜᶜ
function calculate_budget(model, sediment, tracers)
    budget = 0.0
    for tracer in tracers
        budget += sum(getproperty(model.tracers, tracer)[1, 1, 1:model.grid.Nz] .* [Vᶜᶜᶜ(1, 1, k, model.grid) for k in [1:model.grid.Nz;]])
    end
    if sediment
        for sed_tracer in (:Nᵣ, :Nᵣᵣ, :Nᵣₑ)
            budget += getproperty(model.auxiliary_fields, sed_tracer)[1, 1, 1] * (model.grid.Δxᶜᵃᵃ*model.grid.Δyᵃᶜᵃ)
        end
    end
    return budget
end

function calculate_budget(results::OceanBioME.Plot.model_results, grid, bgc_model, sediment)
    budget = zeros(length(results.t))
    for tracer in getproperty(budget_tracers, bgc_model)
        ind = findfirst(results.tracers.=="$tracer")
        budget += sum(results.results[ind, 1, 1, :, :]  .* [Vᶜᶜᶜ(1, 1, k, grid) for k in [1:grid.Nz;]], dims=1)[1, :]
    end
    if sediment
        for sed_tracer in (:Nᵣ, :Nᵣᵣ, :Nᵣₑ)
            ind = findfirst(results.tracers.=="$sed_tracer")
            budget += results.results[ind, 1, 1, 1, :] * (grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ)
        end
    end
    return budget
end

function calculate_C_budget(results::OceanBioME.Plot.model_results, grid, bgc_model, sediment)
    budget = zeros(length(results.t))
    tracers = getproperty(c_budget_tracers, bgc_model)
    for tracer in keys(tracers)
        ind = findfirst(results.tracers.=="$tracer")
        if !isnothing(ind) #incase carbonate model is turned off
            budget += sum(results.results[ind, 1, 1, :, :]  .* [Vᶜᶜᶜ(1, 1, k, grid) for k in [1:grid.Nz;]], dims=1)[1, :] * getproperty(tracers, tracer)
        end
    end

    if sediment
        sed_rd_sym = (:Rdᵣ, :Rdᵣᵣ, :Rd_red)
        for sed_tracer in (:Nᵣ, :Nᵣᵣ, :Nᵣₑ)
            ind = findfirst(results.tracers.=="$sed_tracer")
            budget += results.results[ind, 1, 1, 1, :] * (grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ) * getproperty(OceanBioME.Boundaries.sediment, sed_rd_sym[j])
        end
    end

    airseaparams=merge(Boundaries.defaults.airseaflux, (gas=:CO₂, T=t_function, S=s_function))
    DIC_ind = findfirst(results.tracers.=="DIC")
    ALK_ind = findfirst(results.tracers.=="ALK")
    airsea = zeros(length(results.t))
    for (i, t) in enumerate(results.t)
        budget[i] += Boundaries.airseaflux(.5, .5, t, results.results[DIC_ind, 1, 1, end, i], results.results[ALK_ind, 1, 1, end, i], airseaparams) * (grid.Δxᶜᵃᵃ*grid.Δyᵃᶜᵃ)#if the example does not have a 1x1m x-y grid need to update this to be multipled by area
    end
    return budget
end
end#module