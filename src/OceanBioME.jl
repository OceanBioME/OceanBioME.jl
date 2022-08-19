module OceanBioME

export LOBSTER, NPZ, Light, Boundaries, Particles, Setup, BoxModel, SLatissima

include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Utils/Plot.jl")
include("Models/Biogeochemistry/LOBSTER.jl")
include("Models/Biogeochemistry/NPZ.jl")
include("Models/Macroalgae/SLatissima.jl")
include("Utils/Setup.jl")
include("Utils/BoxModel.jl")

end
