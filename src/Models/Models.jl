module Models

export Sediments

export NPZD, 
       NutrientPhytoplanktonZooplanktonDetritus, 
       LOBSTER,
       PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude

export SLatissima

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
include("Individuals/SLatissima.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")

using .Sediments
using .LOBSTERModel
using .NPZDModel
using .PISCESModel
using .SLatissimaModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
