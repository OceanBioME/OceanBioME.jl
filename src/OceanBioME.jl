"""
A fast and flexible modelling environment for modelling the coupled interactions
between ocean biogeochemistry, carbonate chemistry, and physics.
"""
module OceanBioME

# Biogeochemistry models and useful things
export Biogeochemistry


export NutrientsPlanktonDetritus
export NPZD, LOBSTER, ImplicitBiology, PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

export Nutrients, NitrateAmmonia
export CarbonateSystem, ExplicitCalciumCarbonate
export Abiotic, ImplicitProductivity, PhytoZoo
export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate
export Oxygen

export SugarKelp, SugarKelpParticles, GiantKelp

export InstantRemineralisationSediment, SimpleMultiGSediment
export DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

# Macroalgae models
export BiogeochemicalParticles, SugarKelp, SugarKelpParticles

# Box model
export BoxModel, BoxModelGrid, SpeedyOutput, load_output

# Particles
export Particles

# Light models
export TwoBandPhotosyntheticallyActiveRadiation, 
       PrescribedPhotosyntheticallyActiveRadiation,
       MultiBandPhotosyntheticallyActiveRadiation,
       PrescribedAttenuationPAR,
       PARFromShortwave

# airsea flux
export GasExchange, CarbonDioxideGasExchangeBoundaryCondition, OxygenGasExchangeBoundaryCondition, GasExchangeBoundaryCondition

# carbon chemistry
export CarbonChemistry

# sediment
export Sediments, FlatSediment

# Utilities
export column_advection_timescale, sinking_advection_timescale, Budget

# Positivity preservation utilities
export ClipNegativeTracers, ScaleNegativeTracers, IgnoreNegativeTracerValues

# Oceananigans extensions
export ColumnField, isacolumn

using Oceananigans.Architectures: architecture, device, CPU
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry, AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField
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

struct ContinuousBiogeochemistry{B, L, S, P, M, N} <: AbstractContinuousFormBiogeochemistry
    underlying_biogeochemistry :: B
             light_attenuation :: L
                      sediment :: S
                     particles :: P
                     modifiers :: M
              negative_tracers :: N
end

struct DiscreteBiogeochemistry{B, L, S, P, M, N} <: AbstractBiogeochemistry
    underlying_biogeochemistry :: B
             light_attenuation :: L
                      sediment :: S
                     particles :: P
                     modifiers :: M
              negative_tracers :: N
end

const CompleteBiogeochemistry = Union{<:ContinuousBiogeochemistry, <:DiscreteBiogeochemistry}

@inline resolve_negative_tracers(negative_tracers, underlying_biogeochemistry) = negative_tracers

"""
    Biogeochemistry(underlying_biogeochemistry;
                    light_attenuation = nothing,
                    sediment = nothing,
                    particles = nothing,
                    modifiers = nothing,
                    negative_tracers = nothing)

Construct a biogeochemical model based on `underlying_biogeochemistry` which may have
a `light_attenuation` model, `sediment`, `particles`, `modifiers`, and a `negative_tracers` treatment.

Keyword Arguments
=================

- `light_attenuation`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated
- `negative_tracers`: treatment for negative tracer values; for example [`ClipNegativeTracers`](@ref), [`ScaleNegativeTracers`](@ref), or [`IgnoreNegativeTracerValues`](@ref)
"""
Biogeochemistry(underlying_biogeochemistry;
                light_attenuation = nothing,
                sediment = nothing,
                particles = nothing,
                modifiers = nothing,
                negative_tracers = nothing) = 
    DiscreteBiogeochemistry(underlying_biogeochemistry,
                            light_attenuation,
                            sediment,
                            particles,
                            modifiers,
                            resolve_negative_tracers(negative_tracers, underlying_biogeochemistry))

Biogeochemistry(underlying_biogeochemistry::AbstractContinuousFormBiogeochemistry;
                light_attenuation = nothing,
                sediment = nothing,
                particles = nothing,
                modifiers = nothing,
                negative_tracers = nothing) = 
    ContinuousBiogeochemistry(underlying_biogeochemistry,
                              light_attenuation,
                              sediment,
                              particles,
                              modifiers,
                              resolve_negative_tracers(negative_tracers, underlying_biogeochemistry))

required_biogeochemical_tracers(bgc::CompleteBiogeochemistry) = 
    required_biogeochemical_tracers(bgc.underlying_biogeochemistry)

required_biogeochemical_auxiliary_fields(bgc::CompleteBiogeochemistry) = 
    required_biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry)

biogeochemical_drift_velocity(bgc::CompleteBiogeochemistry, val_name) = 
    biogeochemical_drift_velocity(bgc.underlying_biogeochemistry, val_name)

biogeochemical_auxiliary_fields(bgc::CompleteBiogeochemistry) = 
    merge(biogeochemical_auxiliary_fields(bgc.underlying_biogeochemistry),
          biogeochemical_auxiliary_fields(bgc.light_attenuation))

@inline chlorophyll(bgc::CompleteBiogeochemistry, model) =
    chlorophyll(bgc.negative_tracers, bgc.underlying_biogeochemistry, model)
@inline chlorophyll(::Nothing, model) = ZeroField()

@inline adapt_structure(to, bgc::ContinuousBiogeochemistry) = 
    ContinuousBiogeochemistry(adapt(to, bgc.underlying_biogeochemistry),
                              adapt(to, bgc.light_attenuation),
                              nothing,
                              nothing,
                              nothing,
                              adapt(to, bgc.negative_tracers))

@inline adapt_structure(to, bgc::DiscreteBiogeochemistry) = 
    DiscreteBiogeochemistry(adapt(to, bgc.underlying_biogeochemistry),
                            adapt(to, bgc.light_attenuation),
                            nothing,
                            nothing,
                            nothing,
                            adapt(to, bgc.negative_tracers))

function update_tendencies!(bgc::CompleteBiogeochemistry, model)
    update_tendencies!(bgc, bgc.sediment, model)
    update_tendencies!(bgc, bgc.particles, model)
    update_tendencies!(bgc, bgc.modifiers, model)
end

update_tendencies!(bgc, modifier, model) = nothing
update_tendencies!(bgc, modifiers::Tuple, model) = [update_tendencies!(bgc, modifier, model) for modifier in modifiers]

@inline (bgc::ContinuousBiogeochemistry)(args...) =
    evaluate_continuous_biogeochemistry(bgc.negative_tracers, bgc.underlying_biogeochemistry, args...)

@inline (bgc::DiscreteBiogeochemistry)(i, j, k, grid, val_name, clock, fields) =
    evaluate_discrete_biogeochemistry(bgc.negative_tracers,
                                      bgc.underlying_biogeochemistry,
                                      i, j, k, grid, val_name, clock, fields,
                                      biogeochemical_auxiliary_fields(bgc))

@inline (bgc::ContinuousBiogeochemistry{nothing})(i, j, k, grid, args...) = zero(grid)
@inline (bgc::DiscreteBiogeochemistry{nothing})(i, j, k, grid, args...) = zero(grid)

function update_biogeochemical_state!(bgc::CompleteBiogeochemistry, model)
    # TODO: change the order of arguments here since they should definitly be the other way around
    update_biogeochemical_state!(model, bgc.negative_tracers)
    update_biogeochemical_state!(model, bgc.modifiers)
    update_biogeochemical_state!(model, bgc.light_attenuation)
    update_biogeochemical_state!(model, bgc.underlying_biogeochemistry)
    update_biogeochemical_state!(model, bgc.sediment)
end

update_biogeochemical_state!(model, modifiers::Tuple) = [update_biogeochemical_state!(model, modifier) for modifier in modifiers]

BoxModelGrid(FT = Float64; arch = CPU(), kwargs...) = RectilinearGrid(arch, FT; topology = (Flat, Flat, Flat), size = (), kwargs...)

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
@inline redfield(i, j, k, val_tracer_name, bgc::CompleteBiogeochemistry, tracers) = redfield(i, j, k, val_tracer_name, bgc.underlying_biogeochemistry, tracers)
@inline redfield(val_tracer_name, bgc::CompleteBiogeochemistry) = redfield(val_tracer_name, bgc.underlying_biogeochemistry)
@inline redfield(val_tracer_name, bgc::CompleteBiogeochemistry, tracers) = redfield(val_tracer_name, bgc.underlying_biogeochemistry, tracers)
@inline redfield(val_tracer_name, bgc, tracers) = redfield(val_tracer_name, bgc) 

"""
    conserved_tracers(model::UnderlyingBiogeochemicalModel, args...; kwargs...)

Returns the names of tracers which together are conserved in `model`
"""
conserved_tracers(model::CompleteBiogeochemistry, args...; kwargs...) = conserved_tracers(model.underlying_biogeochemistry, args...; kwargs...)

summary(bgc::CompleteBiogeochemistry) = string("Biogeochemical model based on $(summary(bgc.underlying_biogeochemistry))")
show(io::IO, model::CompleteBiogeochemistry) =
       print(io, summary(model.underlying_biogeochemistry), " \n",
                " Light attenuation: ", summary(model.light_attenuation), "\n",
                " Sediment: ", summary(model.sediment), "\n",
                " Particles: ", summary(model.particles), "\n",
                " Modifiers: ", modifier_summary(model.modifiers), "\n",
                " Negative tracers: ", modifier_summary(model.negative_tracers))

modifier_summary(modifier) = summary(modifier)
modifier_summary(modifiers::Tuple) = "\n"*join([" "*ifelse(n==length(modifiers), "└─", "├─")*summary(modifier) for (n, modifier) in enumerate(modifiers)], "\n")

include("Utils/Utils.jl")
include("Light/Light.jl")
include("Particles/Particles.jl")
include("Sediments/Sediments.jl")
include("BoxModel/boxmodel.jl")
include("Models/Models.jl")

using .Light
using .Particles
using .BoxModels
using .Models

end #module
