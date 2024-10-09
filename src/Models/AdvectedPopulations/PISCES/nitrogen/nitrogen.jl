module Nitrogen

export NitrateAmmonia

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: anoxia_factor, PISCES
using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: oxic_remineralisation, anoxic_remineralisation
using OceanBioME.Models.PISCESModel.Phytoplankton: uptake, nitrogen_availability_limitation, base_production_rate
using OceanBioME.Models.PISCESModel.Zooplankton: upper_trophic_respiration, inorganic_excretion

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

include("nitrate_ammonia.jl")

end # module