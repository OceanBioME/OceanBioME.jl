module Models

export InstantRemineralisationSediment, SimpleMultiGSediment

export NutrientsPlanktonDetritus

export NPZD, LOBSTER, ImplicitBiology, MITgcmDIC, PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

export N, PO₄, Si, Fe
export Nutrients, NitrateAmmonia, SimpleIron
export CarbonateSystem, ExplicitCalciumCarbonate
export Abiotic, ImplicitProductivity, PhytoZoo
export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate
export Oxygen

export SugarKelp, SugarKelpParticles, GiantKelp

export DissolvingMineral, DissolvingParticles, StokesSinking

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
include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl") # NPD's ExplicitCalciumCarbonate uses CarbonChemistry, so load it first
include("AdvectedPopulations/NutrientsPlanktonDetritus/NutrientsPlanktonDetritus.jl")
include("Individuals/SugarKelp/SugarKelp.jl")
include("Individuals/DissolvingParticles/DissolvingParticles.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")

using .SedimentModels
using .NutrientsPlanktonDetritusModels
using .SugarKelpModel
using .DissolvingParticlesModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
