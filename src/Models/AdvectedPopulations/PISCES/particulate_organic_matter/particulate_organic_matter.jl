module ParticulateOrganicMatter

export TwoCompartementCarbonIronParticles

using OceanBioME.Models.PISCESModel: degredation, aggregation, free_iron, PISCES
using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: aggregation_of_colloidal_iron
using OceanBioME.Models.PISCESModel.Phytoplankton: dissolved_exudate, NanoAndDiatoms
using OceanBioME.Models.PISCESModel.Zooplankton: 
    organic_excretion, upper_trophic_excretion, grazing, MicroAndMeso, upper_trophic_fecal_production,
    upper_trophic_fecal_iron_production, calcite_loss

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_drift_velocity
import OceanBioME.Models.PISCESModel: degredation, aggregation
import OceanBioME.Models.PISCESModel.Zooplankton: edible_flux_rate, edible_iron_flux_rate

include("two_size_class.jl")

end # module