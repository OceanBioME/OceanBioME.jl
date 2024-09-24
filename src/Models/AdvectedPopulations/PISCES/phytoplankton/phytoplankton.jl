module Phytoplankton

export NanoAndDiatoms, MixedMondoPhytoplankton, MixedMondoNanoAndDiatoms

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.Zooplankton: grazing

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: mortality

include("nano_and_diatoms.jl")
include("mixed_mondo.jl")
include("growth_rate.jl")
include("nutrient_limitation.jl")
include("mixed_mono_nano_diatoms.jl")

end # module