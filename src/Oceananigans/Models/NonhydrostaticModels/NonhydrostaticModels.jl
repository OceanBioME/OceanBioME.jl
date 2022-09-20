module NonhydrostaticModels

using KernelAbstractions: @index, @kernel, Event, MultiEvent
using KernelAbstractions.Extras.LoopInfo: @unroll

using Oceananigans.Utils: launch!
using Oceananigans.Grids
using Oceananigans.Solvers
using Oceananigans.Forcings: zeroforcing
using Oceananigans: NonhydrostaticModel

import Oceananigans: prognostic_fields

@inline forced_auxiliary_fields(model::NonhydrostaticModel) = Tuple(intersect(keys(model.auxiliary_fields), Tuple(func==zeroforcing ? nothing : c for (c, func) in pairs(model.forcing))))
@inline prognostic_fields(model::NonhydrostaticModel) = merge(model.velocities, model.tracers, NamedTuple{forced_auxiliary_fields(model)}(Tuple(model.auxiliary_fields[c] for c in forced_auxiliary_fields(model))))

include("nonhydrostatic_model.jl")

include("nonhydrostatic_tendency_kernel_functions.jl")
include("calculate_nonhydrostatic_tendencies.jl")
include("set_nonhydrostatic_model.jl")
include("show_nonhydrostatic_model.jl")

end
