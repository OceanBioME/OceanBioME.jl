module BGC

export LOBSTER, NPZ, Light, Boundaries, Particles, Setup, BoxModel

include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Utils/Plot.jl")
include("Models/LOBSTER.jl")
include("Models/NPZ.jl")
include("Utils/Setup.jl")
include("Utils/BoxModel.jl")
end