module Models

export InstantRemineralisationSediment, SimpleMultiGSediment, BurialDenitrificationSediment, BurialDenitrification

export NutrientsPlanktonDetritus

export NPZD, LOBSTER, ImplicitBiology, MARBL, MARBL_Cocco, PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude, PISCESModel

export N, PO₄, Si, Fe
export Nutrients, NitrateAmmonia, ComplexedIron
export CarbonateSystem, ExplicitCalciumCarbonate, ImplicitExplicitCalcite
export Abiotic, ImplicitProductivity, PhytoZoo, ManyPhytoZoo, Autotroph, Heterotroph, Grazing
export MARBL_small_phyto, MARBL_diatoms, MARBL_diazotrophs, MARBL_autotrophs, MARBL_zooplankton,
       MARBL_cocco_small_phyto, MARBL_cocco_diatoms, MARBL_cocco_diazotrophs, MARBL_coccolithophores,
       MARBL_cocco_autotrophs, MARBL_cocco_zooplankton, MARBL_cocco_plankton
export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate, RefractoryDissolvedParticulateCNP, DissolvedRefractoryImplicitCNP, Ballast
export Oxygen, RedoxOxygen

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

include("seawater_density.jl")
include("CarbonChemistry/CarbonChemistry.jl") # NPD's ExplicitCalciumCarbonate uses CarbonChemistry, so load it first
include("AdvectedPopulations/NutrientsPlanktonDetritus/NutrientsPlanktonDetritus.jl")
include("Sediments/Sediments.jl")  # uses the redox helpers, so must follow the NPD models
include("Individuals/SugarKelp/SugarKelp.jl")
include("GasExchange/GasExchange.jl")
include("AdvectedPopulations/PISCES/PISCES.jl")

using .SedimentModels
using .NutrientsPlanktonDetritusModels
using .SugarKelpModel
using .PISCESModel
using .CarbonChemistryModel
using .GasExchangeModel

end # module
