using OceanBioME.Boundaries.Sediments: update_sediment_tendencies!
using OceanBioME.Particles: update_particles_tendencies!

import Oceananigans.Biogeochemistry: update_tendencies!

include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SLatissima.jl")

@inline function update_tendencies!(bgc::ContinuousFormBiogeochemistry, model)
    update_sediment_tendencies!(bgc, bgc.sediment_model, model)
    update_particles_tendencies!(bgc, bgc.particles, model)
    return nothing
end
