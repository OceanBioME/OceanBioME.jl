module OceanBioME

#Biogeochemistry models
export LOBSTER, NPZ, PISCES
#Macroalgae models
export SLatissima
#Setups
export Setup, BoxModel
#Particles
export Particles
#Utilities
export Light, Boundaries, update_timestep!, Budget, Sediments
#Oceananigans extensions
export ColumnField, isacolumn

@inline get_local_value(i, j, k, C) = size(C)[3] == 1 ? C[i, j, 1] : C[i, j, k] #for getting 2D field values

#Overload some Oceananigans functions
include("Oceananigans/Oceananigans.jl")

include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Models/Biogeochemistry/LOBSTER.jl")
include("Models/Biogeochemistry/NPZ.jl")
include("Models/Macroalgae/SLatissima.jl")
include("Utils/Setup.jl")
include("Utils/BoxModel.jl")
include("Utils/Timestep.jl")
include("Utils/Budget.jl")

end
