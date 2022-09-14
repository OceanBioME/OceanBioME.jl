module OceanBioME

#Biogeochemistry models
export LOBSTER, NPZ
#Macroalgae models
export SLatissima
#Setups
export Setup, BoxModel
#Particles
export Particles
#Utilities
export Light, Boundaries, update_timestep!, Budget
#Oceananigans extensions
export ColumnField, isacolumn

#Overload some Oceananigans functions
include("Oceananigans/Oceananigans.jl")

include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Utils/Plot.jl")
include("Models/Biogeochemistry/LOBSTER.jl")
include("Models/Biogeochemistry/NPZ.jl")
include("Models/Macroalgae/SLatissima.jl")
include("Utils/Setup.jl")
include("Utils/BoxModel.jl")
include("Utils/Timestep.jl")
include("Utils/Budget.jl")

end
