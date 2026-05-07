module Models

export InstantRemineralisationSediment, SimpleMultiGSediment

export NPZD, 
       NutrientPhytoplanktonZooplanktonDetritus, 
       LOBSTER, CarbonateSystem, Oxygen, NitrateAmmoniaIron, VariableRedfieldDetritus, Detritus, Nutrient,
       TwoParticleAndDissolved, NitrateAmmonia,
       PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

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

export SimpleCaCO3Precipitation, CaCO3Precipitation

include("Sediments/Sediments.jl")
include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")
include("AdvectedPopulations/CaCO3Precipitation/CaCO3Precipitation.jl")

using .SedimentModels
using .LOBSTERModel
using .NPZDModel
using .SugarKelpModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel
using .CaCO3PrecipitationModel

end # module
