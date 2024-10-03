module Particles 

using Adapt

using KernelAbstractions: @kernel, @index

using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel

using OceanBioME: Biogeochemistry

using Oceananigans.Architectures: architecture, on_architecture

import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Fields: set!
import Oceananigans.Models.LagrangianParticleTracking: update_lagrangian_particle_properties!, step_lagrangian_particles!
import Oceananigans.OutputWriters: fetch_output
import Base: length, size, show, summary
import Adapt: adapt_structure

struct BiogeochemicalParticles{N, B, A, F, T, I, S, X, Y, Z}
      biogeochemistry :: B
            advection :: A
               fields :: F 
          timestepper :: T
  field_interpolation :: I
         scalefactors :: S
                    x :: X # to maintain compatability with lagrangian particle advection from Oceananigans
                    y :: Y
                    z :: Z
end

function required_particle_fields end
function required_tracers end
function coupled_tracers end

@inline required_particle_fields(p::BiogeochemicalParticles) = required_particle_fields(p.biogeochemistry)
@inline required_tracers(p::BiogeochemicalParticles) = required_tracers(p.biogeochemistry)
@inline coupled_tracers(p::BiogeochemicalParticles) = coupled_tracers(p.biogeochemistry)

include("advection.jl")
include("tracer_interpolation.jl")
include("tendencies.jl")
include("time_stepping.jl")
include("update_tracer_tendencies.jl")

function BiogeochemicalParticles(number; 
                                 grid,
                                 biogeochemistry::B,
                                 advection::A = LagrangianAdvection(),
                                 timestepper = ForwardEuler,
                                 field_interpolation::I = NearestPoint(),
                                 scalefactors = ones(number)) where {B, A, I}

    arch = architecture(grid)

    particle_fields = required_particle_fields(biogeochemistry)

    fields = NamedTuple{particle_fields}(ntuple(n->on_architecture(arch, zeros(number)), Val(length(particle_fields))))

    x = on_architecture(arch, zeros(number))
    y = on_architecture(arch, zeros(number))
    z = on_architecture(arch, zeros(number))

    timestepper = timestepper(particle_fields, number, arch)

    scalefactors = on_architecture(arch, scalefactors)

    F = typeof(fields)
    T = typeof(timestepper)
    S = typeof(scalefactors)
    X = typeof(x)
    Y = typeof(y)
    Z = typeof(z)

    return BiogeochemicalParticles{number, B, A, F, T, I, S, X, Y, Z}(biogeochemistry, 
                                                                      advection, 
                                                                      fields, 
                                                                      timestepper, 
                                                                      field_interpolation, 
                                                                      scalefactors,
                                                                      x, y, z)
end

Adapt.adapt_structure(to, p::BiogeochemicalParticles) = 
    BiogeochemicalParticles(adapt(to, p.biogeochemistry),
                            nothing,
                            nothing,
                            nothing,
                            adapt(to, p.field_interpolation),
                            nothing,
                            adapt(to, p.x),
                            adapt(to, p.y),
                            adapt(to, p.z))

# Type piracy...oops
@inline step_lagrangian_particles!(::Nothing, 
    model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, 
                               <:Biogeochemistry{<:Any, <:Any, <:Any, <:BiogeochemicalParticles}}, 
    Δt) = update_lagrangian_particle_properties!(model, model.biogeochemistry, Δt)

@inline step_lagrangian_particles!(::Nothing,
    model::HydrostaticFreeSurfaceModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, 
                                       <:Biogeochemistry{<:Any, <:Any, <:Any, <:BiogeochemicalParticles}, 
                                       <:Any, <:Any, <:Any, <:Any, <:Any,}, 
    Δt) = update_lagrangian_particle_properties!(model, model.biogeochemistry, Δt)

@inline update_lagrangian_particle_properties!(model, bgc::Biogeochemistry{<:Any, <:Any, <:Any, <:BiogeochemicalParticles}, Δt) = 
    update_lagrangian_particle_properties!(bgc.particles, model, bgc, Δt)

@inline function update_lagrangian_particle_properties!(particles::BiogeochemicalParticles, model, bgc, Δt)
    advect_particles!(particles.advection, particles, model, Δt)
    time_step_particle_fields!(particles.timestepper, particles, model, Δt)
end

size(particles::BiogeochemicalParticles) = (length(particles), )
length(::BiogeochemicalParticles{N}) where N = N

# todo, fix below
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

# TODO: implement set!
end #module
