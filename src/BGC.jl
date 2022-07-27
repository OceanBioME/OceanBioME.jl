module BGC

export LOBSTER, NPZ, Light, Boundaries, Particles, Setup, BoxModel

using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using Roots
using Oceananigans.Architectures: device

include("Boundaries.jl")
include("Light.jl")
include("Particles.jl")
include("Plot.jl")
include("LOBSTER.jl")
include("NPZ.jl")
include("Setup.jl")
include("BoxModel.jl")
end