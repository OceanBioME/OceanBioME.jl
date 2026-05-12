module NutrientsPlanktonDetritusModels

export LOBSTER, NPZD,
       CarbonateSystem, 
       Oxygen, 
       NitrateAmmoniaIron, 
       VariableRedfieldDetritus, 
       TwoParticleAndDissolved, 
       NitrateAmmonia,
       Nutrient,
       Detritus,
       NutrientsPlanktonDetritus,
       PhytoZoo

using Oceananigans.Units

using OceanBioME: setup_velocity_fields, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR

import OceanBioME: conserved_tracers, chlorophyll

import Oceananigans.Biogeochemistry: AbstractBiogeochemistry, 
                                     required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

struct NutrientsPlanktonDetritus{NUT, PLA, DET, CAR, OXY} <: AbstractBiogeochemistry
         nutrients :: NUT 
          plankton :: PLA 
          detritus :: DET 
  carbonate_system :: CAR 
            oxygen :: OXY
end

function NutrientsPlanktonDetritus(grid; 
                                   nutrients = nothing,
                                   plankton = nothing,
                                   detritus = nothing,
                                   carbonate_system = nothing,
                                   oxygen = nothing,
                                   light_attenuation = nothing,
                                   sediment = nothing,
                                   scale_negatives = false,
                                   invalid_fill_value = NaN,
                                   particles = nothing,
                                   modifiers = nothing)

    underlying_biogeochemistry = NutrientsPlanktonDetritus(nutrients, plankton, detritus, carbonate_system, oxygen)

    if scale_negatives
        scaler = ScaleNegativeTracers(underlying_biogeochemistry, grid; invalid_fill_value)
        if isnothing(modifiers)
            modifiers = scaler
        elseif modifiers isa Tuple
            modifiers = (modifiers..., scaler)
        else
            modifiers = (modifiers, scaler)
        end
    end
    
    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end

# The possible tracer combinations are:
required_biogeochemical_tracers(bnd::NutrientsPlanktonDetritus) = 
    (required_biogeochemical_tracers(bnd.nutrients)...,
     required_biogeochemical_tracers(bnd.plankton)...,
     required_biogeochemical_tracers(bnd.detritus)...,
     required_biogeochemical_tracers(bnd.carbonate_system)...,
     required_biogeochemical_tracers(bnd.oxygen)...)

# oceananigans defines the fallbacks for ::Nothing or ::AbstractBiogeochemistry but not anything else
required_biogeochemical_tracers(anything_else) = ()
required_biogeochemical_auxiliary_fields(anything_else) = ()

required_biogeochemical_auxiliary_fields(bnd::NutrientsPlanktonDetritus) = 
    (required_biogeochemical_auxiliary_fields(bnd.nutrients)...,
     required_biogeochemical_auxiliary_fields(bnd.plankton)...,
     required_biogeochemical_auxiliary_fields(bnd.detritus)...,
     required_biogeochemical_auxiliary_fields(bnd.carbonate_system)...,
     required_biogeochemical_auxiliary_fields(bnd.oxygen)...)

# fallback - not surer we want this?
@inline (::NutrientsPlanktonDetritus)(i, j, k, grid, val_name, clock, fields, auxiliary_fields) = zero(grid)

include("nutrients.jl")
include("detritus.jl")
include("plankton.jl")
include("carbonate_system.jl")
include("oxygen.jl")
include("show.jl")
include("coupling_utils.jl")
include("adapt_methods.jl")
include("constructors.jl")

end # module