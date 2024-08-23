"""
`CarbonChemistryModel` to solve chemical equilibrium parameterisations
"""
module CarbonChemistryModel

export CarbonChemistry

import Base: show, summary

include("equilibrium_constants.jl")
include("carbon_chemistry.jl")
include("calcite_concentration.jl")

end