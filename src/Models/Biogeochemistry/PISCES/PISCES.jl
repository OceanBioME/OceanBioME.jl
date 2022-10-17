"
PISCES-v2 (Pelagic Interactions Scheme for Carbon and Ecosystem Studies volume 2) is a biogeochemicalmodel  which  simulates  the  lower  trophic  levels  of  marineecosystems  (phytoplankton,  microzooplankton  and  meso-zooplankton) and the biogeochemical cycles of carbon and ofthe main nutrients (P, N, Fe, and Si)Details of the model can be found in the following references: 

References
    (1) O. Aumont1, C. Ethé, A. Tagliabue, L. Bopp, and M. Gehlen, 2015. PISCES-v2: an ocean biogeochemical model for carbon andecosystem studies. Geoscientific Model Development, 8.

Notes
   - All phytoplankton andzooplankton biomasses are in carbon units (mol C L−1) ex-cept  for  the  silicon,  chlorophyll  and  iron  content  of  phy-toplankton,  which  are  respectively  in  Si,  Chl  and  Fe  units(mol Si L−1,  g Chl L−1,  and  mol Fe L−1,  respectively)
   - Convention in PISCES is for z, zₘₓₗ and zₑᵤ to be positive numbers where as in Oceananigans it is conventially a negative number so care must be taken
   - There are several instances throughout where there is some y(x, ...)/x where y(x, ...) = 0 when x = 0, this results in 0/0 = NaN, so following NEMOs solution to this they are modified to y(x, ...)/(x + eps(0.0)) which gives so there is no zero division (eps(0.0) is the machine representation of zero = 5.0e-324)
"
module PISCES
using Oceananigans: Center
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Fields: interpolate, fractional_z_index
using Oceananigans.Operators: ∂zᶜᶜᶜ

include("phytoplankton.jl")
include("zooplankton.jl")
include("doc.jl")
include("two_compartement_POC.jl") #at some point can make this configurable
include("bacteria.jl")
include("nutrients.jl")
include("carbonates.jl")
include("oxygen.jl")

include("par.jl")

@inline get_local_value(i, j, k, C) = size(C)[3] == 1 ? C[i, j, 1] : C[i, j, k] #for getting 2D field values
@inline get_sh(z, zₘₓₗ, shₘₓₗ, shₛᵤ) = ifelse(-zₘₓₗ<z, shₘₓₗ, shₛᵤ)

tracers=(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂)
required_fields = (:PAR¹, :PAR², :PAR³, :T, :S, :zₑᵤ, :zₘₓₗ, :Si̅, :D_dust) #PAR bands, temperature, euphotic layer depth (2D field), mixing layer depth (2D field), maximum Si concentration in a year, dust deposition rate
requried_parameters = (ϕ = "Latitude", L_day = "Function(t) returning length of day ∈ [0, 1]") #latitude (only checked if < 0 to precision not important, except near the equator)

forcing_functions = (P=P_forcing, D=D_forcing, Chlᴾ=Chlᴾ_forcing, Chlᴰ=Chlᴰ_forcing, Feᴾ=Feᴾ_forcing, Feᴰ=Feᴰ_forcing, Siᴰ=Siᴰ_forcing, Z=Z_forcing, M=M_forcing, DOC=DOC_forcing, POC=POC_forcing, GOC=GOC_forcing, Feᴾᴼ=Feᴾᴼ_forcing, Feᴳᴼ=Feᴳᴼ_forcing, Siᴾ=Siᴾ_forcing, NO₃=NO₃_forcing, NH₄=NH₄_forcing, PO₄=PO₄_forcing, Fe=Fe_forcing, Si=Si_forcing, CaCO₃=CaCO₃_forcing, DIC=DIC_forcing, Alk=Alk_forcing, O₂=O₂_forcing)
discrete_forcing = (P=false, D=false, Chlᴾ=false, Chlᴰ=false, Feᴾ=false, Feᴰ=false, Siᴰ=false, Z=false, M=false, DOC=true, POC=true, GOC=true, Feᴾᴼ=true, Feᴳᴼ=true, Siᴾ=true, NO₃=true, NH₄=true, PO₄=true, Fe=true, Si=false, CaCO₃=true, DIC=true, Alk=true, O₂=true)

#parameters
include("parameters.jl")


end # module

#TODO: add Kriest POM model