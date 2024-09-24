module DissolvedOrganicMatter

export DissolvedOrganicCarbon

using OceanBioME.Models.PISCESModel: degredation, aggregation, PISCES
using OceanBioME.Models.PISCESModel.Phytoplankton: dissolved_exudate
using OceanBioME.Models.PISCESModel.Zooplankton: organic_excretion, upper_trophic_excretion

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: degredation, aggregation

include("dissolved_organic_carbon.jl")

end # module