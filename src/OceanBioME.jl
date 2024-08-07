"""
A fast and flexible modelling environment for modelling the coupled interactions
between ocean biogeochemistry, carbonate chemistry, and physics.
"""
module OceanBioME

# Biogeochemistry models
export Biogeochemistry, LOBSTER, NutrientPhytoplanktonZooplanktonDetritus, NPZD, redfield

# Macroalgae models
export SLatissima

# Box model
export BoxModel, BoxModelGrid, SpeedyOutput, load_output

# Particles
export Particles

# Light models
export TwoBandPhotosyntheticallyActiveRadiation, PrescribedPhotosyntheticallyActiveRadiation

# airsea flux
export GasExchange, CarbonDioxideGasExchangeBoundaryCondition, OxygenGasExchangeBoundaryCondition, GasExchangeBoundaryCondition

# carbon chemistry
export CarbonChemistry

# sediment
export Sediments, FlatSediment

# Utilities
export column_advection_timescale, sinking_advection_timescale, Budget

# Positivity preservation utilities
export ScaleNegativeTracers, ZeroNegativeTracers

# Oceananigans extensions
export ColumnField, isacolumn

using Oceananigans.Architectures: architecture, device, CPU
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Grids: RectilinearGrid, Flat

using Adapt
using KernelAbstractions: synchronize

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!,
                                     biogeochemical_drift_velocity,
				                     biogeochemical_auxiliary_fields,
                                     update_tendencies!,
                                     biogeochemical_transition

import Adapt: adapt_structure
import Base: show, summary

struct Biogeochemistry{B, L, S, P, M} <: AbstractContinuousFormBiogeochemistry
    underlying_biogeochemistry :: B
             light_attenuation :: L
                      sediment :: S
                     particles :: P
                     modifiers :: M
    
    Biogeochemistry(underlying_biogeochemistry::B,
                    light_attenuation::L,
                    sediment::S,
                    particles::P,
                    modifiers::M) where {B, L, S, P, M} = 
        new{B, L, S, P, M}(underlying_biogeochemistry,
                           light_attenuation,
                           sediment,
                           particles,
                           modifiers)
end

"""
    Biogeochemistry(underlying_biogeochemistry;
                    light_attenuation = nothing,
                    sediment = nothing,
                    particles = nothing,
                    modifiers = nothing)

Construct a biogeochemical model based on `underlying_biogeochemistry` which may have
a `light_attenuation` model, `sediment`, `particles`, and `modifiers`.

Keyword Arguments
=================

- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated
"""
Biogeochemistry(underlying_biogeochemistry;
                light_attenuation = nothing,
                sediment = nothing,
                particles = nothing,
                modifiers = nothing) = 
    Biogeochemistry(underlying_biogeochemistry,
                    light_attenuation,
                    sediment,
                    particles,
                    modifiers)

required_biogeochemical_tracers(bgc::Biogeochemistry) = required_biogeochemical_tracers(bgc.underlying_biogeochemistry)

required_biogeochemical_auxiliary_fields(bgc::Biogeochemistry) = required_biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry)

biogeochemical_drift_velocity(bgc::Biogeochemistry, val_name) = biogeochemical_drift_velocity(bgc.underlying_biogeochemistry, val_name)

biogeochemical_auxiliary_fields(bgc::Biogeochemistry) = merge(biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry),
                                                              biogeochemical_auxiliary_fields(bgc.light_attenuation))

@inline adapt_structure(to, bgc::Biogeochemistry) = adapt(to, bgc.underlying_biogeochemistry)

function update_tendencies!(bgc::Biogeochemistry, model)
    update_tendencies!(bgc, bgc.sediment, model)
    update_tendencies!(bgc, bgc.particles, model)
    update_tendencies!(bgc, bgc.modifiers, model)
    synchronize(device(architecture(model)))
end

update_tendencies!(bgc, modifier, model) = nothing
update_tendencies!(bgc, modifiers::Tuple, model) = [update_tendencies!(bgc, modifier, model) for modifier in modifiers]

# do we still need this for CPU kernels???
@inline biogeochemical_transition(i, j, k, grid, bgc::Biogeochemistry, val_tracer_name, clock, fields) =
    biogeochemical_transition(i, j, k, grid, bgc.underlying_biogeochemistry, val_tracer_name, clock, fields) 

@inline (bgc::Biogeochemistry)(args...) = bgc.underlying_biogeochemistry(args...)

function update_biogeochemical_state!(bgc::Biogeochemistry, model)
    update_biogeochemical_state!(model, bgc.modifiers)
    synchronize(device(architecture(model)))
    update_biogeochemical_state!(model, bgc.light_attenuation)
end

update_biogeochemical_state!(model, modifiers::Tuple) = [update_biogeochemical_state!(model, modifier) for modifier in modifiers]

BoxModelGrid(; arch = CPU(), kwargs...) = RectilinearGrid(arch; topology = (Flat, Flat, Flat), size = (), kwargs...)

@inline maximum_sinking_velocity(bgc) = 0.0

"""
    redfield(i, j, k, val_tracer_name, bgc, tracers)

Returns the redfield ratio of `tracer_name` from `bgc` at `i`, `j`, `k`.
"""
@inline redfield(i, j, k, val_tracer_name, bgc, tracers) = redfield(val_tracer_name, bgc, tracers)

"""
    redfield(val_tracer_name, bgc)
    redfield(val_tracer_name, bgc, tracers)

Returns the redfield ratio of `tracer_name` from `bgc` when it is constant across the domain.
"""
@inline redfield(val_tracer_name, bgc) = NaN

# fallbacks
@inline redfield(i, j, k, val_tracer_name, bgc::Biogeochemistry, tracers) = redfield(i, j, k, val_tracer_name, bgc.underlying_biogeochemistry, tracers)
@inline redfield(val_tracer_name, bgc::Biogeochemistry) = redfield(val_tracer_name, bgc.underlying_biogeochemistry)
@inline redfield(val_tracer_name, bgc::Biogeochemistry, tracers) = redfield(val_tracer_name, bgc.underlying_biogeochemistry, tracers)
@inline redfield(val_tracer_name, bgc, tracers) = redfield(val_tracer_name, bgc) 

"""
    conserved_tracers(model::UnderlyingBiogeochemicalModel)

Returns the names of tracers which together are conserved in `model`
"""
conserved_tracers(model::Biogeochemistry) = conserved_tracers(model.underlying_biogeochemistry)

summary(bgc::Biogeochemistry) = string("Biogeochemical model based on $(summary(bgc.underlying_biogeochemistry))")
show(io::IO, model::Biogeochemistry) =
       print(io, summary(model.underlying_biogeochemistry), " \n",
                " Light attenuation: ", summary(model.light_attenuation), "\n",
                " Sediment: ", summary(model.sediment), "\n",
                " Particles: ", summary(model.particles), "\n",
                " Modifiers: ", summary(model.modifiers))

summary(modifiers::Tuple) = tuple([summary(modifier) for modifier in modifiers])

include("Utils/Utils.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("BoxModel/boxmodel.jl")
include("Models/Models.jl")

using .Light
using .BoxModels
using .Models

end #module
