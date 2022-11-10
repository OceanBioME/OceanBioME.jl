"
Light attenuation by chlorophyll as described by [Karleskind2011](@cite) (implimented as twoBand) and [Morel1988](@cite).

These should be used by setting up a callback like:
`simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)), TimeStepCallsite())`
"
module Light

using Oceananigans
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device

include("2band.jl")
include("morel.jl")

end