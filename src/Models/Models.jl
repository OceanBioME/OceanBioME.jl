module Models

export InstantRemineralisationSediment, SimpleMultiGSediment

export NutrientsPlanktonDetritus

export NPZD, LOBSTER, ImplicitBiology, PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

export N, PO₄, Si, Fe
export Nutrients, NitrateAmmonia
export CarbonateSystem, ExplicitCalciumCarbonate
export Abiotic, ImplicitProductivity, PhytoZoo
export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate
export Oxygen

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
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl") # NPD's ExplicitCalciumCarbonate uses CarbonChemistry, so load it first
include("AdvectedPopulations/NutrientsPlanktonDetritus/NutrientsPlanktonDetritus.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")
include("AdvectedPopulations/CaCO3Precipitation/CaCO3Precipitation.jl")

using .SedimentModels
using .NutrientsPlanktonDetritusModels
using .SugarKelpModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel
using .CaCO3PrecipitationModel

end # module
