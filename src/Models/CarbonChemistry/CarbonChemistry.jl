"""
`CarbonChemistryModel` to solve chemical equilibrium parameterisations
"""
module CarbonChemistryModel

export CarbonChemistry

include("equilibrium_constants.jl")
include("carbonate_chemistry.jl")

end