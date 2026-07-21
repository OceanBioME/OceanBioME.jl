module Models

export InstantRemineralisationSediment, SimpleMultiGSediment

export NutrientsPlanktonDetritus

export NPZD, LOBSTER, ImplicitBiology, PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

export N, PO₄, Si, Fe
export Nutrients, NitrateAmmonia
export CarbonateSystem
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

include("Sediments/Sediments.jl")
include("AdvectedPopulations/NutrientsPlanktonDetritus/NutrientsPlanktonDetritus.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")

using .SedimentModels
using .NutrientsPlanktonDetritusModels
using .SugarKelpModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
