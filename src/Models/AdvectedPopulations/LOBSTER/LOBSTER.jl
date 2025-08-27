module LOBSTERModel

export LOBSTER

using Oceananigans.Units

using OceanBioME: setup_velocity_fields, Biogeochemistry
using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR

import Oceananigans.Biogeochemistry: AbstractBiogeochemistry, 
                                     required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

# default to "standard" LOBSTER
@kwdef struct LOBSTER{NUT, BIO, DET, CAR, OXY} <: AbstractBiogeochemistry
         nutrients :: NUT = NitrateAmmonia() # NitrateAmmoniaIron()
           biology :: BIO = PhytoZoo() # 
          detritus :: DET = TwoParticleAndDissolved() # VariableRedfield()
  carbonate_system :: CAR = nothing # CarbonateSystem()
            oxygen :: OXY = nothing
end

function LOBSTER(grid; 

                 nutrients = NitrateAmmonia(),
                 biology = PhytoZoo(),
                 detritus = TwoParticleAndDissolved(grid),
                 carbonate_system = nothing,
                 oxygen = nothing,

                 surface_photosynthetically_active_radiation = default_surface_PAR,

                 light_attenuation =
                       TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                                  surface_PAR = surface_photosynthetically_active_radiation),
                   
                 sediment = nothing,

                 scale_negatives = false,
                 invalid_fill_value = NaN,

                 particles = nothing,
                 modifiers = nothing)

    if !isnothing(sediment) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    lobster = LOBSTER(; biology, nutrients, detritus, carbonate_system, oxygen)

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
    
    return Biogeochemistry(lobster;
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end

# The possible tracer combinations are:
required_biogeochemical_tracers(lobster::LOBSTER) = 
    (required_biogeochemical_tracers(lobster.nutrients)...,
     required_biogeochemical_tracers(lobster.biology)...,
     required_biogeochemical_tracers(lobster.detritus)...,
     required_biogeochemical_tracers(lobster.carbonate_system)...,
     required_biogeochemical_tracers(lobster.oxygen)...)

required_biogeochemical_tracers(::Nothing) = ()

# we should never need this maybe
#@inline (::LOBSTER)(args...) = 0 # fallback for `Nothing` models

include("nutrients.jl")
include("biology.jl")
include("detritus.jl")
include("carbonate_system.jl")
include("oxygen.jl")

end # module