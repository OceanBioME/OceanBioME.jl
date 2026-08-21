# MARBLS nPmZ plankton model where phytoplankton pools have carbon,
# chlorophyll, phosphorus, iron, and optionally silicon and calcite
# and zooplankton pools have carbon only
module MARBLPlankton # maybe rename

using Oceananigans: defaults
using OceanBioME.Models.NutrientsPlanktonDetritusModels.PlanktonModels: AbstractPlankton

include("autotroph.jl")
include("heterotroph.jl")
include("primatives.jl")
include("many_phyto_zoo.jl")
include("infrastructure.jl")
include("manifest_many_phyto_zoo.jl")

end # module