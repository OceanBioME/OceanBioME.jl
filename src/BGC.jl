module BGC

export Lobster, Light, AirSeaFlux

using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using Roots
using Oceananigans.Architectures: device

mutable struct BGCModel{T, F, B}
    tracers :: T
    forcing :: F
    boundary_conditions :: B
end

include("AirSeaFlux.jl")
include("Light.jl")
include("Lobster.jl")

end