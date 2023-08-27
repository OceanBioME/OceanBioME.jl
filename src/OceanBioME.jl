"""
A fast and flexible modelling environment for modelling the coupled interactions
between ocean biogeochemistry, carbonate chemistry, and physics.
"""
module OceanBioME

# Biogeochemistry models
export LOBSTER, NutrientPhytoplanktonZooplanktonDetritus, PISCES, NPZD, redfield

# Macroalgae models
export SLatissima

# Box model
export BoxModel, BoxModelGrid, SaveBoxModel, run!, set!

# Particles
export Particles

# Light models
export TwoBandPhotosyntheticallyActiveRadiation

# Boundaries
export Boundaries, Sediments, GasExchange, FlatSediment

# Utilities
export column_advection_timescale, column_diffusion_timescale, sinking_advection_timescale, Budget

# Positivity preservaiton utilities
export zero_negative_tracers!, error_on_neg!, warn_on_neg!, ScaleNegativeTracers, remove_NaN_tendencies!

# Oceananigans extensions
export ColumnField, isacolumn

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Adapt

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!,
                                     biogeochemical_drift_velocity,
				                     biogeochemical_auxiliary_fields,
                                     update_tendencies!

import Adapt: adapt_structure
import Base: show, summary

struct Biogeochemistry{B, L, S, P, I} <: AbstractContinuousFormBiogeochemistry
    underlying_biogeochemistry :: B
             light_attenuation :: L
                      sediment :: S
                     particles :: P
                        inputs :: I
    
    Biogeochemistry(underlying_biogeochemistry::B,
                    light_attenuation::L,
                    sediment::S,
                    particles::P,
                    inputs::I) where {B, L, S, P, I} = 
        new{B, L, S, P, I}(underlying_biogeochemistry,
                           light_attenuation,
                           sediment,
                           particles,
                           inputs)
end

"""
    Biogeochemistry(underlying_biogeochemistry;
                    light_attenuation = nothing,
                    sediment = nothing,
                    particles = nothing,
                    inputs = nothing)

Construct a biogeochemical model based on `underlying_biogeochemistry` which may have
a `light_attenuation` model, `sediment`, `particles`, and `inputs`.

Keyword Arguments
=================

- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `particles`: slot for `BiogeochemicalParticles`
- `inputs`: slot for nutrient inputs such as rivers (work in progress)
"""
Biogeochemistry(underlying_biogeochemistry;
                light_attenuation = nothing,
                sediment = nothing,
                particles = nothing,
                inputs = nothing) = 
    Biogeochemistry(underlying_biogeochemistry,
                    light_attenuation,
                    sediment,
                    particles,
                    inputs)

required_biogeochemical_tracers(bgc::Biogeochemistry) = required_biogeochemical_tracers(bgc.underlying_biogeochemistry)

required_biogeochemical_auxiliary_fields(bgc::Biogeochemistry) = required_biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry)

biogeochemical_drift_velocity(bgc::Biogeochemistry, val_name) = biogeochemical_drift_velocity(bgc.underlying_biogeochemistry, val_name)

biogeochemical_auxiliary_fields(bgc::Biogeochemistry) = merge(biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry),
                                                              biogeochemical_auxiliary_fields(bgc.light_attenuation))

adapt_structure(to, bgc::Biogeochemistry) = Biogeochemistry(adapt(to, bgc.underlying_biogeochemistry),
                                                            adapt(to, bgc.light_attenuation),
                                                            adapt(to, bgc.sediment),
                                                            adapt(to, bgc.particles),
                                                            adapt(to, bgc.inputs))

function update_tendencies!(bgc::Biogeochemistry, model)
    update_tendencies!(bgc, bgc.sediment, model)
    update_tendencies!(bgc, bgc.particles, model)
    update_tendencies!(bgc, bgc.inputs, model)
end

@inline (bgc::Biogeochemistry)(args...) = bgc.underlying_biogeochemistry(args...)

function update_biogeochemical_state!(bgc::Biogeochemistry, model)
    update_biogeochemical_state!(model, bgc.light_attenuation)
end

update_tendencies!(bgc, ::Nothing, model) = nothing

abstract type UnderlyingBiogeochemicalModel end

@inline (::UnderlyingBiogeochemicalModel)(val_tracer_name, x, y, z, t, fields...) = zero(x)

struct BoxModelGrid end

@inline maximum_sinking_velocity(bgc) = 0.0

"""
    redfield(i, j, k, val_tracer_name, bgc, tracers)

Returns the redfield ratio of `tracer_name` from `bgc` at `i`, `j`, `k`.
"""
@inline redfield(i, j, k, val_tracer_name, bgc, tracers) = NaN

"""
    redfield(val_tracer_name, bgc, tracers)

Returns the redfield ratio of `tracer_name` from `bgc` when it is constant across the domain.
"""
@inline redfield(val_tracer_name, bgc) = NaN

@inline redfield(i, j, k, val_tracer_name, bgc::Biogeochemistry, tracers) = redfield(i, j, k, val_tracer_name, bgc.underlying_biogeochemistry, tracers)
@inline redfield(val_tracer_name, bgc::Biogeochemistry) = redfield(val_tracer_name, bgc.underlying_biogeochemistry)

summary(bgc::Biogeochemistry) = string("Biogeochemical model based on $(summary(bgc.underlying_biogeochemistry))")
show(io::IO, model::Biogeochemistry) =
       print(io, show(model.underlying_biogeochemistry), " \n",
                " Light attenuation: ", summary(model.light_attenuation), "\n",
                " Sediment: ", summary(model.sediment), "\n",
                " Particles: ", summary(model.particles))

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
import .SLatissimaModel.SLatissima

end #module
