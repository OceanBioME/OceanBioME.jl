module DissolvedOrganicMatter

export DissolvedOrganicCarbon

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: 
    degredation, aggregation, PISCES, free_iron, anoxia_factor
using OceanBioME.Models.PISCESModel.Phytoplankton: dissolved_exudate
using OceanBioME.Models.PISCESModel.Zooplankton: 
    organic_excretion, upper_trophic_excretion, bacteria_concentration, bacteria_activity

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: degredation, aggregation

include("dissolved_organic_carbon.jl")

end # module