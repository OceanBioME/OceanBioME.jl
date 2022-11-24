module OceanBioME

# Biogeochemistry models
export LOBSTER, NutrientPhytoplanktonZooplanktonDetritus, PISCES
# Macroalgae models
export SLatissima
# Setups
export Setup, BoxModel
# Particles
export Particles
# Light models
export TwoBandPhotosyntheticallyActiveRatiation, update_PAR!
# Utilities
export Boundaries, update_timestep!, Budget, Sediments, GasExchange
# Positivity preservaiton utilities
export zero_negative_tracers!, error_on_neg!, warn_on_neg!, scale_negative_tracers!
# Oceananigans extensions
export ColumnField, isacolumn

@inline get_local_value(i, j, k, C) = size(C)[3] == 1 ? C[i, j, 1] : C[i, j, k] #for getting 2D field values

include("Utils/Utils.jl")
include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Models/Models.jl")

using .Boundaries
using .Light

end
