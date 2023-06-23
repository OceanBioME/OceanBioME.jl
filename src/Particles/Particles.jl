module Particles 

using OceanBioME: ContinuousFormBiogeochemistry
using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel

import Oceananigans.Models.LagrangianParticleTracking: update_lagrangian_particle_properties!, step_lagrangian_particles!
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.OutputWriters: fetch_output
import Base: length, size, show, summary

abstract type BiogeochemicalParticles end

# TODO: add model.particles passing

@inline step_lagrangian_particles!(::Nothing, 
    model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, 
                               <:ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}}, 
    Δt) = update_lagrangian_particle_properties!(model, model.biogeochemistry, Δt)

@inline step_lagrangian_particles!(::Nothing,
    model::HydrostaticFreeSurfaceModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, 
                                       <:ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, 
                                       <:Any, <:Any, <:Any, <:Any, <:Any,}, 
    Δt) = update_lagrangian_particle_properties!(model, model.biogeochemistry, Δt)

@inline update_lagrangian_particle_properties!(model, bgc::ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, Δt) = 
    update_lagrangian_particle_properties!(bgc.particles, model, bgc, Δt)

update_lagrangian_particle_properties!(::BiogeochemicalParticles, model, bgc, Δt) = nothing
update_tendencies!(bgc::ContinuousFormBiogeochemistry{<:Any, <:Any, <:BiogeochemicalParticles}, model) = update_tendencies!(bgc, bgc.particles, model)
update_tendencies!(bgc, ::BiogeochemicalParticles, model) = nothing

size(particles::BiogeochemicalParticles) = size(particles.x)
length(particles::BiogeochemicalParticles) = length(particles.x)

Base.summary(particles::BiogeochemicalParticles) =
    string(length(particles), " BiogeochemicalParticles with eltype ", nameof(eltype(particles)),
           " and properties ", propertynames(particles))

function Base.show(io::IO, particles::BiogeochemicalParticles)
    Tparticle = nameof(eltype(particles))
    properties = propertynames(particles)
    Nparticles = length(particles)

    print(io, Nparticles, " BiogeochemicalParticles with eltype ", Tparticle, ":", "\n",
        "└── ", length(properties), " properties: ", properties, "\n")
end

# User may want to overload this to not output parameters over and over again
function fetch_output(particles::BiogeochemicalParticles, model)
    names = propertynames(particles)
    return NamedTuple{names}([getproperty(particles, name) for name in names])
end

include("tracer_tendencies.jl")
end#module