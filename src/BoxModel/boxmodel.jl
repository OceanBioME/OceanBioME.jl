"
Integrates the biogeochemical models in a closed box
"
module BoxModels

export BoxModel, run!, set!, SaveBoxModel

using Oceananigans.Biogeochemistry: 
        AbstractContinuousFormBiogeochemistry, 
        required_biogeochemical_tracers, 
        required_biogeochemical_auxiliary_fields

using Oceananigans: Clock, prettytime
using Oceananigans.TimeSteppers: tick!

using OceanBioME: BoxModelGrid
using StructArrays, JLD2

import Oceananigans.Simulations: run!
import Oceananigans: set!

@inline no_func(args...) = 0.0

mutable struct BoxModel{B, V, FT, F, TS, C}
    biogeochemistry :: B
    values :: V
    stop_time :: FT
    forcing :: F
    timestepper :: TS
    Δt :: FT
    clock :: C
end

function BoxModel(;biogeochemistry::B,
                   stop_time::FT = 0.0,
                   forcing = NamedTuple(),
                   timestepper::TS = RungeKutta3TimeStepper((required_biogeochemical_tracers(biogeochemistry)..., required_biogeochemical_auxiliary_fields(biogeochemistry)...)),
                   Δt::FT = 1.0,
                   clock::C = Clock(0.0, 0, 1)) where {B, FT, TS, C}

    variables = (required_biogeochemical_tracers(biogeochemistry)..., required_biogeochemical_auxiliary_fields(biogeochemistry)...)
    values = StructArray([NamedTuple{variables}(zeros(FT, length(variables)))])
    forcing = NamedTuple{variables}([variable in keys(forcing) ? forcing[variable] : no_func for variable in variables])

    V = typeof(values)
    F = typeof(forcing)

    return BoxModel{B, V, FT, F, TS, C}(biogeochemistry, values, stop_time, forcing, timestepper, Δt, clock)
end

function calculate_tendencies!(model)
    # TODO: add sinking out of box
    Gⁿ = model.timestepper.Gⁿ

    for tracer in required_biogeochemical_tracers(model.biogeochemistry)
        @inbounds getproperty(Gⁿ, tracer) .= model.biogeochemistry(Val(tracer), 0.0, 0.0, 0.0, model.clock.time, model.values[1]...) + model.forcing[tracer](model.clock.time, model.values[1]...)
    end

    for variable in required_biogeochemical_auxiliary_fields(model.biogeochemistry)
        @inbounds getproperty(Gⁿ, variable) .= model.forcing[variable](model.clock.time, model.values[1]...)
    end
end

# a user could overload this like `update_boxmodel_state!(model::BoxModel{<:SomeBiogeochemicalModel, <:Any, <:Any, <:Any, <:Any, <:Any})`
# this could be used e.g. to update a PAR field
@inline update_boxmodel_state!(model::BoxModel) = nothing

# should abstract out to simulation like Oceananians to add e.g. callbacks
function run!(model::BoxModel; feedback_interval = 1000, save_interval = Inf, save = nothing)
    itter = 0
    while model.clock.time < model.stop_time
        itter % feedback_interval == 0 && @info "Reached $(prettytime(model.clock.time))"
        time_step!(model, model.Δt)
        itter += 1

        itter % save_interval == 0 && save(model)
    end

    if !isnothing(save)
        close(save.file)
    end
end

function set!(model::BoxModel; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ keys(model.values[1])
            # hacky with StructArrays
            getproperty(model.values, fldname) .= value
        else
            throw(ArgumentError("name $fldname not found in model."))
        end
    end
end

struct SaveBoxModel{FP, F}
    filepath :: FP
    file :: F

    function SaveBoxModel(filepath::FP) where FP
        file = jldopen(filepath, "w+")
        F = typeof(file)
        return new{FP, F}(filepath, file)
    end
end

(save::SaveBoxModel)(model) = save.file["values/$(model.clock.time)"] = model.values[1]

include("timesteppers.jl")

end # module