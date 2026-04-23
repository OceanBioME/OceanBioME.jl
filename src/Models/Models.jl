module Models

export InstantRemineralisationSediment, SimpleMultiGSediment

export NPZD, LOBSTER, CarbonateSystem, Oxygen, NitrateAmmoniaIron, VariableRedfieldDetritus, Detritus, Nutrient,
       TwoParticleAndDissolved, NitrateAmmonia,
       PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel,
       BiologyNutrientDetritus, PhytoZoo

export SugarKelp, SugarKelpParticles, GiantKelp

export CarbonChemistry

export GasExchange, 
       CarbonDioxideGasExchangeBoundaryCondition, 
       OxygenGasExchangeBoundaryCondition, 
       GasExchangeBoundaryCondition,
       ScaledGasTransferVelocity,
       SchmidtScaledTransferVelocity,
       CarbonDioxidePolynomialSchmidtNumber,
       OxygenPolynomialSchmidtNumber

include("Sediments/Sediments.jl")
include("AdvectedPopulations/BiologyNutrientDetritus/BiologyNutrientDetritus.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")

using .SedimentModels
using .BiologyNutrientDetritusModels
using .SugarKelpModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
