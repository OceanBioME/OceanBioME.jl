module NutrientsPlanktonDetritusModels

export NutrientsPlanktonDetritus

export LOBSTER, NPZD, ImplicitBiology

export Nutrients, N, PO₄, Si, Fe, NitrateAmmonia
export CarbonateSystem
export Abiotic, ImplicitProductivity, PhytoZoo
export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate
export Oxygen

using Adapt
using Oceananigans.Grids: AbstractGrid

import Adapt: adapt_structure
import Base: show, summary

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

# Consider all biogeochemical models to abstract to "Nutrients" (inorganics which can limit)
# growth, "Plankton" which is all the planktonic biology (typically phyto and zoo, but 
# could be mixotrophs), and "detritus" which is all inorganic waste. You then also have 
# "inorganic carbon" and "oxygen" which should just plug in to the others (they are essentially
# one way coupled to the rest of the system, kind of).

# We (some what arbitarily) make the abstract flow between these groups be:
# Nutrients -> Plankton: `nutrient_uptake` with `nutrient_limitation` effecting the plankton
# Plankton -> Nutrients: `inorganic_waste` or `inorganic_X_waste`
# Plankton -> Detritus: `solid_waste` and `dissolved_waste` or ...`_X_waste`
# Detritus -> Plankton: `grazing` which can limit plankton
# Detritus -> Nutrients: `inorganic_waste` or `inorganic_X_waste`

# we are going to assume that the default units are nitrogen, which is maybe a wrong assumtion to some

include("nutrients_plankton_detritus.jl")
include("assumptions.jl")

include("Nutrients/Nutrients.jl")
include("InorganicCarbon/InorganicCarbon.jl")
include("Detritus/Detritus.jl")
include("Plankton/Plankton.jl")
include("Oxygen/Oxygen.jl")

using .NutrientsModels
using .InorganicCarbonModels
using .PlanktonModels
using .DetritusModels
using .OxygenModels

include("constructor.jl")
include("utils.jl")

end # module