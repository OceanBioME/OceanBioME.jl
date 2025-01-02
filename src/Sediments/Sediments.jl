module Sediments


#### !!! This should all be a set of boundary conditions !!!

using KernelAbstractions: @kernel, @index

using Oceananigans: Clock
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: TimeStepper

import Oceananigans.Models: fields, prognostic_fields

import Oceananigans.Biogeochemistry: update_biogeochemical_state!, 
                                     update_tendencies!

struct BiogeochemicalSediment{BC, TS, CL, GR, SF, TF, BI} <: AbstractModel{TS}
    biogeochemistry :: BC
        timestepper :: TS
              clock :: CL
               grid :: GR
             fields :: SF
     tracked_fields :: TF
     bottom_indices :: BI
end

@inline fields(s::BiogeochemicalSediment) = s.fields
@inline prognostic_fields(s::BiogeochemicalSediment) = fields(s) # this can have different methods if auxiliary fields are required

@inline required_sediment_fields(s::BiogeochemicalSediment) = required_sediment_fields(s.biogeochemistry)
@inline required_tracers(s::BiogeochemicalSediment) = required_tracers(s.biogeochemistry)
@inline sinking_fluxs(s::BiogeochemicalSediment) = sinking_fluxs(s.biogeochemistry)
@inline coupled_tracers(s::BiogeochemicalSediment) = coupled_tracers(s.biogeochemistry)

abstract type AbstractSedimentBiogeochemistry end
abstract type AbstractContinuousFormSedimentBiogeochemistry <: AbstractSedimentBiogeochemistry end

const ACFSBGC = AbstractContinuousFormSedimentBiogeochemistry

include("timesteppers.jl")
include("bottom_indices.jl")
include("tracked_fields.jl")
include("update_state.jl")
include("compute_tendencies.jl")
include("tracer_coupling.jl")

function BiogeochemicalSediment(grid, biogeochemistry;
                                clock = Clock(time = zero(grid)),
                                timestepper = :QuasiAdamsBashforth2)

    bottom_indices = calculate_bottom_indices(grid)

    sediment_field_names = required_sediment_fields(biogeochemistry)
    tracked_field_names = (required_tracers(biogeochemistry)..., sinking_fluxs(biogeochemistry)...)

    prognostic_fields = NamedTuple{sediment_field_names}(map(n -> Field{Center, Center, Nothing}(grid), 1:length(sediment_field_names)))
    tracked_fields = NamedTuple{tracked_field_names}(map(n -> Field{Center, Center, Nothing}(grid), 1:length(tracked_field_names)))

    timestepper = TimeStepper(timestepper, grid, prognostic_fields)

    return BiogeochemicalSediment(biogeochemistry, timestepper, clock, grid, prognostic_fields, tracked_fields, bottom_indices)
end

end # module