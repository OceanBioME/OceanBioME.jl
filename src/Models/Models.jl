module Models

export Sediments

export NPZD, 
       NutrientPhytoplanktonZooplanktonDetritus, 
       LOBSTER,
       PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel,
       MARBL

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
include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")
include("AdvectedPopulations/MARBL/MARBL.jl")

using .Sediments
using .LOBSTERModel
using .NPZDModel
using .SugarKelpModel
using .PISCESModel
using .MARBLModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
