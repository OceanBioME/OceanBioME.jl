module SedimentModels

export InstantRemineralisation, InstantRemineralisationSediment,
       SimpleMultiG, SimpleMultiGSediment

using Adapt

using OceanBioME.Sediments: AbstractContinuousFormSedimentBiogeochemistry,
                            BiogeochemicalSediment

import OceanBioME.Sediments: required_sediment_fields,
                             required_tracers,
                             sinking_fluxs,
                             coupled_tracers

import Adapt: adapt_structure

include("simple_multi_G.jl")
include("instant_remineralisation.jl")

end # module