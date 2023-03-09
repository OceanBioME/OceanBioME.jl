module Particles 

using OceanBioME: ContinuousFormBiogeochemistry
using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel

import Oceananigans.LagrangianParticleTracking: update_particle_properties!
import Oceananigans.Biogeochemistry: update_tendencies!
import Base: length

abstract type BiogeochemicalParticles end

@inline update_particle_properties!(particles, 
                                    model::NonhydrostaticModel{TS, E, A, G, T, B, R, SD, U, C, Φ, F, V, S, K, BG, P, 
                                                               ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, 
                                                               I, AF}, 
                                    Δt) where {TS, E, A, G, T, B, R, SD, U, C, Φ, F, V, S, K, BG, P, I, AF} =
        update_particle_properties!(particles, model, model.biogeochemistry, Δt)

@inline update_particle_properties!(particles, 
                                    model::HydrostaticFreeSurfaceModel{TS, E, A, S, G, T, V, B, R, F, P, 
                                                                       ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, 
                                                                       U, C, Φ, K, AF}, 
                                    Δt) where {TS, E, A, S, G, T, V, B, R, F, P, U, C, Φ, K, AF} =
        update_particle_properties!(particles, model, model.biogeochemistry, Δt)

@inline update_particle_properties!(particles, ::Any, ::ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, Δt) = 
    error("We currently do not support biogeochemical particles along side other particles.
You may be able to impliment your needs along side as custom dynamics of the 
biogeochemical particles but if not please make an issue and we will try and 
resolve this problem.")

@inline update_particle_properties!(::Nothing, ::Any, ::ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, Δt) = 
    update_particle_properties!(bgc.particles, model, bgc, Δt)

update_particle_properties!(::BiogeochemicalParticles, model, bgc, Δt) = nothing
update_tendencies!(bgc::ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, model) = update_tendencies!(bgc, bgc.particles, model)
update_tendencies!(bgc, ::BiogeochemicalParticles, model) = nothing

@inline length(particles::BiogeochemicalParticles) = length(particles.x)

include("tracer_tendencies.jl")
end#module
