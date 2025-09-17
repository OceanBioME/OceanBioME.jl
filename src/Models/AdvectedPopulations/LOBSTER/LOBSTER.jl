"""
The Lodyc-DAMTP Ocean Biogeochemical Simulation Tools for Ecosystem and
Resources (LOBSTER) model.

LOBSTER is an NPZD type model with multiple nutrients and detritus pools.
In the default configuration there is one phytoplankton and one zooplankton,
nitrate and ammonia, and small, large, and dissovled detritus.

Optionally it can also include iron (which limits phytoplankton growth),
the carbonate system (DIC and alkalinity), oxygen, and variable N:C ratio
detritus where each of dissolved, small, and large particles have both 
nitrogen and carbon compartements.

Required submodels
==================

* Photosynthetically available radiation: PAR (W/m²)

When carbonate system is active:
* Temperature: T (ᵒC)
* Salinity: S (‰)
"""
module LOBSTERModel

export LOBSTER, CarbonateSystem, Oxygen, NitrateAmmoniaIron, VariableRedfieldDetritus

using Oceananigans.Units

using OceanBioME: setup_velocity_fields, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR

import OceanBioME: conserved_tracers, chlorophyll_ratio

import Oceananigans.Biogeochemistry: AbstractBiogeochemistry, 
                                     required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity
                                     

# default to "standard" LOBSTER
struct LOBSTER{NUT, BIO, DET, CAR, OXY} <: AbstractBiogeochemistry
         nutrients :: NUT # = NitrateAmmonia() # NitrateAmmoniaIron()
           biology :: BIO # = PhytoZoo() # 
          detritus :: DET # = TwoParticleAndDissolved() # VariableRedfieldDetritus()
  carbonate_system :: CAR # = nothing # CarbonateSystem()
            oxygen :: OXY # = nothing
end

"""
    LOBSTER(; grid

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

Construct an instance of the [LOBSTER](@ref LOBSTER) biogeochemical model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to configure sinking speeds
- `nutrients`: nutrint component, defaults to `NitrateAmmonia` which has nitrate (`NO₃`) and ammonia (`NH₄`),
currently `NitrateAmmoniaIron` is the alternative which includes iron (`Fe`)
- `biology`: biological (living) component, defaults to `PhytoZoo` which has phytoplankton (`P`) and zooplankton (`Z`)
- `detritus`: non-living organic component, default to `TwoParticlesAndDissolved` which has small and large 
particles (`sPOM` and `bPOM`) and dissolved organic (`DOM`), alternativly can be `VariableRedfieldDetritus` which has 
`sPON`, `bPON`, `DON`, `sPOC`, `bPOC`, and `DOC`
- `carbonate_system`: an optional model component for the inorganic carbon components, defaults to `nothing` but
can be `CarbonateSystem` which has dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`)
- `oxygen`: an optional model component for oxygen, defaults to `nothing` but can be `Oxygen` which has oxygen (`O₂`)
- `surface_photosynthetically_active_radiation`: funciton for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `BiogeochemicalSediment`
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> model = LOBSTER(; grid)
LOBSTER model (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM) 
 Light attenuation: Two-band light attenuation model (Float64)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing

```
"""
function LOBSTER(; grid, 

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

    lobster = LOBSTER(nutrients, biology, detritus, carbonate_system, oxygen)

    if scale_negatives
        scaler = ScaleNegativeTracers(lobster, grid; invalid_fill_value)
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
include("show.jl")
include("coupling_utils.jl")

end # module