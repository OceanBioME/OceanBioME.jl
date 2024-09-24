module Zooplankton

export MicroAndMesoZooplankton, QualityDependantZooplankton, MicroAndMeso

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: anoxia_factor, PISCES, flux_rate

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: mortality

function edible_flux_rate end
function edible_iron_flux_rate end

include("food_quality_dependant.jl")
include("micro_and_meso.jl")
include("defaults.jl")

end # module