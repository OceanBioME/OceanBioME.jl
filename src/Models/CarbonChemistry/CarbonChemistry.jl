"""
`CarbonChemistryModel` to solve chemical equilibrium parameterisations
"""
module CarbonChemistryModel

export CarbonChemistry

import Base: show, summary

const ATM = 101325 # Pa
const GAS_CONSTANT = 8.31446261815324 # J / kg / mol

include("equilibrium_constants.jl")
include("virial_coefficients.jl")
include("carbon_chemistry.jl")
include("calcite_concentration.jl")

end