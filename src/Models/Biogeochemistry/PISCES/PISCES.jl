"
PISCES-v2 (Pelagic Interactions Scheme for Car-bon and Ecosystem Studies volume 2) is a biogeochemicalmodel  which  simulates  the  lower  trophic  levels  of  marineecosystems  (phytoplankton,  microzooplankton  and  meso-zooplankton) and the biogeochemical cycles of carbon and ofthe main nutrients (P, N, Fe, and Si)Details of the model can be found in the following references: 

References
    (1) O. Aumont1, C. Ethé, A. Tagliabue, L. Bopp, and M. Gehlen, 2015. PISCES-v2: an ocean biogeochemical model for carbon andecosystem studies. Geoscientific Model Development, 8.

Notes
    All phytoplankton andzooplankton biomasses are in carbon units (mol C L−1) ex-cept  for  the  silicon,  chlorophyll  and  iron  content  of  phy-toplankton,  which  are  respectively  in  Si,  Chl  and  Fe  units(mol Si L−1,  g Chl L−1,  and  mol Fe L−1,  respectively)
    Convention in PISCES is for z, zₘₓₗ and zₑᵤ to be positive numbers where as in Oceananigans it is conventially a negative number so care must be taken
"
module PISCES
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Fields: interpolate
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

@inline get_local_value(i, j, k, C) = ifelse(length(size(C)) == 3, C[i, j, k], C[i, j, 1]) #for getting 2D field values
@inline get_sh(z, zₘₓₗ, shₘₓₗ, shₛᵤ) = ifelse(-zₘₓₗ<z, shₘₓₗ, shₛᵤ)

tracers=(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂)
required_fields = (:PAR¹, :PAR², :PAR³, :T, :S, :zₑᵤ, :zₘₓₗ, :ϕ) #PAR bands, temperature, euphotic layer depth (2D field), mixing layer depth (2D field), latitude (2D field)

forcing_functions = (P=P_forcing, D=D_forcing, Chlᴾ=Chlᴾ_forcing, Chlᴰ=Chlᴰ_forcing, Feᴾ=Feᴾ_forcing, Feᴰ=Feᴰ_forcing, Siᴰ=Siᴰ_forcing, Z=Z_forcing, M=M_forcing, DOC=DOC_forcing, POC=POC_forcing, GOC=GOC_forcing, Feᴾᴼ=Feᴾᴼ_forcing, Feᴳᴼ=Feᴳᴼ_forcing, Siᴾ=Siᴾ_forcing, NO₃=NO₃_forcing, NH₄=NH₄_forcing, PO₄=PO₄_forcing, Fe=Fe_forcing, Si=Si_forcing, CaCO₃=CaCO₃_forcing, DIC=DIC_forcing, Alk=Alk_forcing, O₂=O₂_forcing)
discrete_forcing = (P=false, D=false, Chlᴾ=false, Chlᴰ=false, Feᴾ=false, Feᴰ=false, Siᴰ=false, Z=false, M=false, DOC=true, POC=true, GOC=true, Feᴾᴼ=true, Feᴳᴼ=true, Siᴾ=true, NO₃=true, NH₄=true, PO₄=false, Fe=false, Si=false, CaCO₃=true, DIC=false, Alk=true, O₂=true)

#parameters
include("parameters.jl")
end # module

#TODO: add Kriest POM model

#= Interesting way to prevent negative tracers used in NEMOs PISCES implimentation

        Handling of the negative concentrations
         ! The biological SMS may generate negative concentrations
         ! Trends are tested at each grid cell. If a negative concentrations 
         ! is created at a grid cell, all the sources and sinks at that grid 
         ! cell are scale to avoid that negative concentration. This approach 
         ! is quite simplistic but it conserves mass.

          xnegtr(:,:,:) = 1.e0
         DO jn = jp_pcs0, jp_pcs1
            DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, jpk)
               IF( ( tr(ji,jj,jk,jn,Kbb) + tr(ji,jj,jk,jn,Krhs) ) < 0.e0 ) THEN
                  ztra             = ABS( tr(ji,jj,jk,jn,Kbb) ) / ( ABS( tr(ji,jj,jk,jn,Krhs) ) + rtrn )
                  xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
               ENDIF
            END_3D
         END DO
         !                                ! where at least 1 tracer concentration becomes negative
         !                                ! 
         ! Concentrations are updated
         DO jn = jp_pcs0, jp_pcs1
           tr(:,:,:,jn,Kbb) = tr(:,:,:,jn,Kbb) + xnegtr(:,:,:) * tr(:,:,:,jn,Krhs)
         END DO

         =#