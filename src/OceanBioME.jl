module OceanBioME

# Biogeochemistry models
export LOBSTER, NutrientPhytoplanktonZooplanktonDetritus, PISCES

# Macroalgae models
export SLatissima

# Box model
export BoxModel, BoxModelGrid, SaveBoxModel, run!, set!

# Particles
export Particles

# Light models
export TwoBandPhotosyntheticallyActiveRatiation, update_PAR!

# Boundaries
export Boundaries, Sediments, GasExchange

# Utilities
export column_advection_timescale, column_diffusion_timescale, sinking_adveciton_timescale, Budget

# Positivity preservaiton utilities
export zero_negative_tracers!, error_on_neg!, warn_on_neg!, ScaleNegativeTracers

# Oceananigans extensions
export ColumnField, isacolumn

@inline get_local_value(i, j, k, C) = size(C)[3] == 1 ? C[i, j, 1] : C[i, j, k] #for getting 2D field values

struct BoxModelGrid end

@inline maximum_sinking_velocity(bgc) = 0.0

include("Utils/Utils.jl")
include("Boundaries/Boundaries.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("BoxModel/boxmodel.jl")
include("Models/Models.jl")

using .Boundaries
using .Light
using .BoxModels
using .LOBSTERModel
using .NPZDModel

end
