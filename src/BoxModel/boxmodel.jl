"""
Integrate biogeochemical models on a single point
"""
module BoxModels

export BoxModel, run!, set!, SpeedyOutput, load_output

using Oceananigans: Clock, prettytime
using Oceananigans.Biogeochemistry: 
        AbstractContinuousFormBiogeochemistry, 
        required_biogeochemical_tracers, 
        required_biogeochemical_auxiliary_fields,
        update_biogeochemical_state!

using Oceananigans.Fields: CenterField
using Oceananigans.TimeSteppers: tick!, TimeStepper

using OceanBioME: BoxModelGrid
using StructArrays, JLD2

import Oceananigans.Simulations: run!
import Oceananigans: set!
import Oceananigans.Fields: regularize_field_boundary_conditions, TracerFields
import Oceananigans.Architectures: architecture
import Oceananigans.Models: default_nan_checker, iteration, AbstractModel, prognostic_fields
import Oceananigans.TimeSteppers: update_state!
import Oceananigans.OutputWriters: default_included_properties

import Base: show, summary

@inline no_func(args...) = 0.0

mutable struct BoxModel{G, B, F, FV, FO, TS, C, PT} <: AbstractModel{TS}
               grid :: G # here so that simualtion can be built
    biogeochemistry :: B
             fields :: F
       field_values :: FV
            forcing :: FO
        timestepper :: TS
              clock :: C
  prescribed_tracers :: PT
end

"""
    BoxModel(; biogeochemistry,
               forcing = NamedTuple(),
               timestepper = :RungeKutta3,
               clock = Clock(; time = 0.0),
               prescribed_tracers::PT = (T = (t) -> 0, ))

Constructs a box model of a `biogeochemistry` model. Once this has been constructed you can set initial condiitons by `set!(model, X=1.0...)`.

Keyword Arguments 
=================

- `biogeochemistry`: (required) an OceanBioME biogeochemical model, most models must be passed a `grid` which can be set to a `BoxModelGrid` for box models
- `forcing`: NamedTuple of additional forcing functions for the biogeochemical tracers to be integrated
- `timestepper`: Timestepper to integrate model
- `clock`: Oceananigans clock to keep track of time
- `prescribed_tracers`: named tuple of tracer names and function (`f(t)`) prescribing tracer values
"""
function BoxModel(; biogeochemistry::B,
                    grid = BoxModelGrid(),
                    forcing = NamedTuple(),
                    timestepper = :RungeKutta3,
                    clock::C = Clock(; time = 0.0),
                    prescribed_tracers::PT = (T = (t) -> 0, )) where {B, C, PT}

    variables = required_biogeochemical_tracers(biogeochemistry)
    fields = NamedTuple{variables}([CenterField(grid) for var in eachindex(variables)])
    field_values = zeros(length(variables)+length(required_biogeochemical_auxiliary_fields(biogeochemistry)))
    forcing = NamedTuple{variables}([variable in keys(forcing) ? forcing[variable] : no_func for variable in variables])

    timestepper = BoxModelTimeStepper(timestepper, grid, keys(fields))

    F = typeof(fields)
    FV = typeof(field_values)
    FO = typeof(forcing)
    TS = typeof(timestepper)

    return BoxModel{typeof(grid), B, F, FV, FO, TS, C, PT}(grid, biogeochemistry, fields, field_values, forcing, timestepper, clock, prescribed_tracers)
end

function update_state!(model::BoxModel, callbacks=[]; compute_tendencies = true)
    t = model.clock.time

    for (name, forcing) in pairs(model.prescribed_tracers)
        if name in keys(model.fields)
            @inbounds model.fields[name][1, 1, 1] = @inbounds forcing(t)
        end
    end

    for callback in callbacks
        callback.callsite isa UpdateStateCallsite && callback(model)
    end

    update_biogeochemical_state!(model.biogeochemistry, model)

    compute_tendencies && 
        compute_tendencies!(model, callbacks)

    return nothing
end

architecture(model::BoxModel) = architecture(model.grid) # this might be the default fallback
default_nan_checker(::BoxModel) = nothing
iteration(model::BoxModel) = model.clock.iteration
prognostic_fields(model::BoxModel) = @inbounds model.fields[required_biogeochemical_tracers(model.biogeochemistry)]


"""
    set!(model::BoxModel; kwargs...)

Set the `values` for a `BoxModel`

Arguments
=========

- `model` - the model to set the arguments for

Keyword Arguments
==================

- variables and value pairs to set
"""
function set!(model::BoxModel; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(model.fields)
            ϕ = getproperty(model.fields, fldname)
        else
            throw(ArgumentError("name $fldname not found in model.fields."))
        end
        set!(ϕ, value)
    end
    return nothing
end

default_included_properties(::BoxModel) = [:grid]

include("timesteppers.jl")
include("output_writer.jl")

summary(::BoxModel)  = string("Biogeochemical box model")
show(io::IO, model::BoxModel{G, B, F, FV, FO, TS, C, PT}) where {G, B, F, FV, FO, TS, C, PT} = 
       print(io, summary(model), "\n",
                "  Biogeochemical model: ", "\n",
                "    └── ", summary(model.biogeochemistry), "\n",
                "  Time-stepper:", "\n", 
                "    └── ", nameof(typeof(model.timestepper)), "\n",
                "  Time:", "\n",
                "    └── $(prettytime(model.clock.time))")

end # module
