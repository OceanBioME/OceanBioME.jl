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
module BiologyNutrientDetritusModels

export LOBSTER, 
       CarbonateSystem, 
       Oxygen, 
       NitrateAmmoniaIron, 
       VariableRedfieldDetritus, 
       TwoParticleAndDissolved, 
       NitrateAmmonia,
       Nutrient,
       Detritus

using Oceananigans.Units

using OceanBioME: setup_velocity_fields, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR

import OceanBioME: conserved_tracers, chlorophyll

import Oceananigans.Biogeochemistry: AbstractBiogeochemistry, 
                                     required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity
                                     

struct BiologyNutrientDetritus{NUT, BIO, DET, CAR, OXY} <: AbstractBiogeochemistry
         nutrients :: NUT 
           biology :: BIO 
          detritus :: DET 
  carbonate_system :: CAR 
            oxygen :: OXY 
end

# The possible tracer combinations are:
required_biogeochemical_tracers(bnd::BiologyNutrientDetritus) = 
    (required_biogeochemical_tracers(bnd.nutrients)...,
     required_biogeochemical_tracers(bnd.biology)...,
     required_biogeochemical_tracers(bnd.detritus)...,
     required_biogeochemical_tracers(bnd.carbonate_system)...,
     required_biogeochemical_tracers(bnd.oxygen)...)

# needed?
required_biogeochemical_tracers(::Nothing) = ()

required_biogeochemical_auxiliary_fields(bnd::BiologyNutrientDetritus) = 
    (required_biogeochemical_auxiliary_fields(bnd.nutrients)...,
     required_biogeochemical_auxiliary_fields(bnd.biology)...,
     required_biogeochemical_auxiliary_fields(bnd.detritus)...,
     required_biogeochemical_auxiliary_fields(bnd.carbonate_system)...,
     required_biogeochemical_auxiliary_fields(bnd.oxygen)...)

# fallback - not surer we want this?
@inline (::BiologyNutrientDetritus)(i, j, k, grid, val_name, clock, fields, auxiliary_fields) = zero(grid)

include("nutrients.jl")
include("detritus.jl")
include("biology.jl")
include("carbonate_system.jl")
include("oxygen.jl")
include("show.jl")
include("coupling_utils.jl")
include("adapt_methods.jl")

end # module