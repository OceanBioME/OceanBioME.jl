"""
`CarbonChemistryModel` to solve chemical equilibrium parameterisations
"""
module CarbonChemistryModel

export CarbonChemistry

import Base: show, summary

include("equilibrium_constants.jl")
include("carbonate_chemistry.jl")

end