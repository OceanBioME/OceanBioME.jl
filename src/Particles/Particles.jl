module Particles 

export BiogeochemicalParticles

using Adapt

using KernelAbstractions: @kernel, @index

using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel
using OceanBioME: DiscreteBiogeochemistry, ContinuousBiogeochemistry

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: AbstractGrid

import Oceananigans.Architectures: architecture
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Fields: set!
import Oceananigans.Models.LagrangianParticleTracking: update_lagrangian_particle_properties!, step_lagrangian_particles!
import Oceananigans.OutputWriters: fetch_output
import Base: length, size, show, summary
import Adapt: adapt_structure

abstract type AbstractBiogeochemicalParticles end

struct BiogeochemicalParticles{N, B, A, F, T, I, S, X, Y, Z} <: AbstractBiogeochemicalParticles
      biogeochemistry :: B
            advection :: A
               fields :: F 
          timestepper :: T
  field_interpolation :: I
         scalefactors :: S
                    x :: X # to maintain compatability with lagrangian particle advection from Oceananigans
                    y :: Y
                    z :: Z

    BiogeochemicalParticles(number, biogeochemistry::B, advection::A,
                                    fields::F, timestepper::T, 
                                    field_interpolation::I, scalefactors::S,
                                    x::X, y::Y, z::Z) where {B, A, F, T, I, S, X, Y, Z} =
        new{number, B, A, F, T, I, S, X, Y, Z}(biogeochemistry, advection, fields, timestepper,
                                               field_interpolation, scalefactors, x, y, z)     
end

architecture(p::BiogeochemicalParticles) = architecture(p.x)

function required_particle_fields end
function required_tracers end
function coupled_tracers end

@inline required_particle_fields(p::BiogeochemicalParticles) = required_particle_fields(p.biogeochemistry)
@inline required_tracers(p::BiogeochemicalParticles) = required_tracers(p.biogeochemistry)
@inline coupled_tracers(p::BiogeochemicalParticles) = coupled_tracers(p.biogeochemistry)

include("atomic_operations.jl")
include("advection.jl")
include("tracer_interpolation.jl")
include("tendencies.jl")
include("time_stepping.jl")
include("update_tracer_tendencies.jl")
include("set.jl")

"""
    BiogeochemicalParticles(number; 
                            grid,
                            biogeochemistry,
                            advection = LagrangianAdvection(),
                            timestepper = ForwardEuler,
                            field_interpolation = NearestPoint(),
                            scalefactors = ones(number))

Creates `number` particles with `biogeochemistry` on `grid`, advected by
`advection` which defaults to `LagrangianAdvection` (i.e. they comove with
the water). The biogeochemistry is stepped by `timestepper` and tracer fields
are interpolated by `field_interpolation`, which defaults to directly reading
the nearest center point and taking up/depositing in the same.

Particles can also have a `scalefactor` which scales their tracer interaction
(e.g. to mimic the particle representing multiple particles).
"""
function BiogeochemicalParticles(number; 
                                 grid::AbstractGrid{FT},
                                 biogeochemistry,
                                 advection = LagrangianAdvection(),
                                 timestepper = ForwardEuler,
                                 field_interpolation = NearestPoint(),
                                 scalefactors = ones(number)) where FT

    arch = architecture(grid)

    particle_fields = required_particle_fields(biogeochemistry)

    fields = NamedTuple{particle_fields}(ntuple(n->on_architecture(arch, zeros(FT, number)), Val(length(particle_fields))))

    x = on_architecture(arch, zeros(FT, number))
    y = on_architecture(arch, zeros(FT, number))
    z = on_architecture(arch, zeros(FT, number))

    timestepper = timestepper(particle_fields, number, arch)

    scalefactors = on_architecture(arch, scalefactors)

    return BiogeochemicalParticles(number, 
                                   biogeochemistry, 
                                   advection, 
                                   fields, 
                                   timestepper, 
                                   field_interpolation, 
                                   scalefactors,
                                   x, y, z)
end

Adapt.adapt_structure(to, p::BiogeochemicalParticles{N}) where N = 
    BiogeochemicalParticles(N, adapt(to, p.biogeochemistry),
                               nothing,
                               adapt(to, p.fields),
                               nothing,
                               adapt(to, p.field_interpolation),
                               adapt(to, p.scalefactors),
                               adapt(to, p.x),
                               adapt(to, p.y),
                               adapt(to, p.z))

const BGC_WITH_PARTICLES = Union{<:DiscreteBiogeochemistry{<:Any, <:Any, <:Any, <:AbstractBiogeochemicalParticles},
                                 <:ContinuousBiogeochemistry{<:Any, <:Any, <:Any, <:AbstractBiogeochemicalParticles}}

const NONHYDROSTATIC_WITH_BGC_PARTICLES = NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any,
                                                              <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any,
                                                              <:Any, <:Any, <:Any, <:Any, <:BGC_WITH_PARTICLES}

const HYDROSTATIC_WITH_BGC_PARTICLES = HydrostaticFreeSurfaceModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any,
                                                                   <:Any, <:Any, <:Any, <:Any, <:Any,
                                                                   <:BGC_WITH_PARTICLES}

const MODEL_WITH_BGC_PARTICLES = Union{<:NONHYDROSTATIC_WITH_BGC_PARTICLES, <:HYDROSTATIC_WITH_BGC_PARTICLES}

@inline step_lagrangian_particles!(::Nothing, model::MODEL_WITH_BGC_PARTICLES, Δt) =
    step_lagrangian_particles!(model.biogeochemistry, model, Δt)

@inline step_lagrangian_particles!(bgc::BGC_WITH_PARTICLES, model, Δt) =
    update_lagrangian_particle_properties!(bgc.particles, model, bgc, Δt)

@inline function update_lagrangian_particle_properties!(particles::BiogeochemicalParticles, model, bgc, Δt)
    advect_particles!(particles.advection, particles, model, Δt)
    time_step_particle_fields!(particles.timestepper, particles, model, Δt)
    update_particle_state!(particles, model, Δt)
end

@inline update_particle_state!(particles, model, Δt) = nothing

size(particles::BiogeochemicalParticles) = (length(particles), )
length(::BiogeochemicalParticles{N}) where N = N

# todo, fix below
Base.summary(particles::BiogeochemicalParticles{N}) where N =
    string(N, " BiogeochemicalParticles with ", summary(particles.biogeochemistry))

function Base.show(io::IO, particles::BiogeochemicalParticles)
    properties = required_particle_fields(particles)
    tracers = coupled_tracers(particles)

    print(io, summary(particles), ":", "\n",
        "├── fields: ", properties, "\n",
        "└── coupled tracers: ", tracers, "\n")
end

# User may want to overload this to not output parameters over and over again
function fetch_output(particles::BiogeochemicalParticles, model)
    names = propertynames(particles.fields)
    return NamedTuple{names}([getproperty(particles.fields, name) for name in names])
end

end #module
