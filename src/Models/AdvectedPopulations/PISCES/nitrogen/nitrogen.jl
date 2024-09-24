module Nitrogen

export NitrateAmmonia

using OceanBioME.Models.PISCESModel: anoxia_factor, PISCES
using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: oxic_remineralisation, anoxic_remineralisaiton
using OceanBioME.Models.PISCESModel.Phytoplankton: uptake, nitrogen_availability_limitation, base_production_rate
using OceanBioME.Models.PISCESModel.Zooplankton: upper_trophic_respiration, inorganic_excretion

include("nitrate_ammonia.jl")

end # module