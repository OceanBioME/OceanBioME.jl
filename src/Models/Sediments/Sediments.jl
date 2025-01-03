module SedimentModels

export InstantRemineralisation, InstantRemineralisationSediment,
       SimpleMultiG, SimpleMultiGSediment

using OceanBioME.Sediments: AbstractContinuousFormSedimentBiogeochemistry,
                            BiogeochemicalSediment

import OceanBioME.Sediments: required_sediment_fields,
                             required_tracers,
                             sinking_fluxs,
                             coupled_tracers


include("simple_multi_G.jl")
include("instant_remineralization.jl")

end # module