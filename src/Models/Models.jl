using OceanBioME.Particles: BiogeochemicalParticles
using OceanBioME.Boundaries: AbstractSediment

import Oceananigans.Biogeochemistry: update_tendencies!

include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SLatissima.jl")

@inline function update_tendencies!(bgc::ContinuousFormBiogeochemistry{<:Any, <:AbstractSediment, <:BiogeochemicalParticles}, model)
    update_tendencies!(bgc, bgc.sediment_model, model)
    update_tendencies!(bgc, bgc.particles, model)
end