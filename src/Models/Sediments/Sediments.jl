module SedimentModels

export InstantRemineralisation, InstantRemineralisationSediment,
       SimpleMultiG, SimpleMultiGSediment,
       BurialDenitrification, BurialDenitrificationSediment

using Adapt
using Oceananigans.Units: day
using Oceananigans.Operators: Δzᶜᶜᶠ
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell

using OceanBioME.Sediments: AbstractContinuousFormSedimentBiogeochemistry,
                            AbstractSedimentBiogeochemistry,
                            BiogeochemicalSediment

using ..NutrientsPlanktonDetritusModels.OxygenModels: suboxic_fraction, aerobic_fraction

import OceanBioME.Sediments: required_sediment_fields,
                             required_tracers,
                             sinking_fluxes,
                             coupled_tracers

import Adapt: adapt_structure

include("simple_multi_G.jl")
include("instant_remineralisation.jl")
include("burial_denitrification.jl")
include("adapt_show_methods.jl")

end # module