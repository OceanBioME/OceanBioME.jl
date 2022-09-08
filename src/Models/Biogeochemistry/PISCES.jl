"
PISCES-v2 (Pelagic Interactions Scheme for Car-bon and Ecosystem Studies volume 2) is a biogeochemicalmodel  which  simulates  the  lower  trophic  levels  of  marineecosystems  (phytoplankton,  microzooplankton  and  meso-zooplankton) and the biogeochemical cycles of carbon and ofthe main nutrients (P, N, Fe, and Si)Details of the model can be found in the following references: 

References
    (1) O. Aumont1, C. Ethé, A. Tagliabue, L. Bopp, and M. Gehlen, 2015. PISCES-v2: an ocean biogeochemical model for carbon andecosystem studies. Geoscientific Model Development, 8.

Notes
    All phytoplankton andzooplankton biomasses are in carbon units (mol C L−1) ex-cept  for  the  silicon,  chlorophyll  and  iron  content  of  phy-toplankton,  which  are  respectively  in  Si,  Chl  and  Fe  units(mol Si L−1,  g Chl L−1,  and  mol Fe L−1,  respectively)
"
#module PISCES
using Oceananigans
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

@inline zₘₓₗ(x, y, params) = @inbounds params.zₘₓₗ[findfirst(params.x.==x), findfirst(params.y.==y)]
@inline zₑᵤ(x, y, params) = @inbounds params.zₑᵤₗ[findfirst(params.x.==x), findfirst(params.y.==y)]
@inline ϕ(y, params::NamedTuple) = ϕ(y, params.ϕ)
@inline ϕ(y, lat::Float64) = lat
@inline ϕ(y, lat::Function) = lat(y) #lat being an interpolation function
#Probably more options for what lat could be

@inline μᴾ(L_day, zₘₓₗ, zₑᵤ, P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, PARᴾ, T, I, params) = μₚ(T, params)*f₁(L_day, params)*f₂(zₘₓₗ, zₑᵤ, I, params)*(1-exp(-getproperty(params.α, I)*(Chlᴾ/P)*PARᴾ/(L_day*(params.μ_ref+b_resp))))*Lₗᵢₘ(P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, I, params)

@inline μₚ(T, I, params) = params.μ⁰ₘₐₓ*getproperty(params.b, I)^T

@inline f₁(L_day, params) = 1.5*L_day/(0.5+L_day)

@inline t_dark(zₘₓₗ, zₑᵤ) = max(0, zₘₓₗ-zₑᵤ)^2/86400
@inline f₂(zₘₓₗ, zₑᵤ, I::Symbol, params) = 1 - t_dark(zₘₓₗ, zₑᵤ)/(t_dark(zₘₓₗ, zₑᵤ)+getproperty(params.t_dark, I))

@inline L_NO₃(NO₃, NH₄, I, params) = getproperty(params.K_NO₃, I)*NH₄/(getproperty(params.K_NO₃, I)*getproperty(params.K_NH₄, I)+getproperty(params.K_NH₄, I)*NO₃+getproperty(params.K_NO₃, I)*NH₄)
@inline L_NH₄(NO₃, NH₄, I, params) = getproperty(params.K_NH₄, I)*NO₃/(getproperty(params.K_NO₃, I)*getproperty(params.K_NH₄, I)+getproperty(params.K_NH₄, I)*NO₃+getproperty(params.K_NO₃, I)*NH₄)
@inline Lₙ(NO₃, NH₄, I, params) = L_NO₃(NO₃, NH₄, I, params) + L_NO₃(NO₃, NH₄, I, params) + L_NH₄(NO₃, NH₄, I, params)
@inline Lₚₒ₄(PO₄, P, I, params) = PO₄/(PO₄+K(P, I, :PO₄, params))

@inline L_Fe(P, Chlᴾ, Feᴾ, NO₃, NH₄, I, params) = min(1, max(0, ((Feᴾ/P) - θᶠᵉₘᵢₙ(P, Chlᴾ, NO₃, NH₄, I, params))/getproperty(params.θᶠᵉₒₚₜ, I)))

@inline Lₗᵢₘ(P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, I, params) = min(Lₚₒ₄(PO₄, P, I, params), 
                                                                            Lₙ(NO₃, NH₄, I, params),
                                                                            L_Fe(p, Chlᴾ, Feᴾ, NO₃, NH₄, I, params),
                                                                            I==:D ? params.Kₛᵢᴰᵐⁱⁿ + 7*params.Si̅^2/(params.Kₛᵢ^2+params.Si̅^2) : Inf)

@inline K(P, I, J, params) = getproperty(getproperty(params.Kᵐⁱⁿ, I), J)*(P₁(P, I, params)+getproperty(params.Sᵣₐₜ, I)*P₂(P, I, params))/(P₁(P, I, params)+P₂(P, I, params))
@inline P₁(P, I, J, params) = min(P, getproperty(params.Pₘₐₓ, I))
@inline P₂(P, I, J, params) = min(0, P - getproperty(params.Pₘₐₓ, I))

P_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :P, params)
D_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :D, params)

@inline Phyto_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, I::Symbol, params) = (
    (1-getproperty(params.δ, I))*μᴾ(params.L_day(t), zₘₓₗ(x, y, params), zₑᵤ(x, y, params), (I==:P ? [P, Chlᴾ, Feᴾ] : [D, Chlᴰ, Feᴰ])..., NO₃, NH₄, PO₄, Fe, I==:P ? PARᴾ : PARᴰ, T, I, params)*P 
    -getproperty(params.m, I)*(I==:P ? P : D)/(params.Kₘ+(I==:P ? P : D))*(I==:P ? P : D)
    -(zₘₓₗ(x, y, params)<z ? params.shₘₓₗ : params.shₛᵤ)*(I==:P ? params.wᴾ : params.wᵖ + params.wₘₐₓᴰ*(1-Lₗᵢₘ(D, Feᴰ, NO₃, NH₄, PO₄, Fe, :D, params)))*(I==:P ? P : D)^2
    -g(P, D, POC, I, :Z, params)*Z-g(P, D, POC, I, :M, params)*M)

@inline F(P, D, POC, J::Symbol, params) = getproperty(params.p, J).P*max(0, P-getproperty(params.Jₜₕᵣ, J).P)+getproperty(params.p, J).D*max(0, D-getproperty(params.Jₜₕᵣ, J).D)+getproperty(params.p, J).POC*max(0, POC-getproperty(params.Jₜₕᵣ, J).POC)
@inline Fₗᵢₘ(P, D, POC, J::Symbol, params) = max(0, F(P, D, POC, J, params)-min(0.5*F(P, D, POC, J, params), getproperty(params.Fₜₕᵣ, J)))
@inline g(P, D, POC, I::Symbol, J::Symbol, params) = getproperty(params.gₘ, J)*(Fₗᵢₘ(P, D, POC, J, params)/F(P, D, POC, J, params))*getproperty(getproperty(params.p, J), I)*max(0, (I == :P ? P : (I == D ? D : POC)) - getproperty(getproperty(params.Zₜₕᵣ, J), I))/(getproperty(params.K_G, J)+getproperty(params.p, J).P*P+getproperty(params.p, J).D*D+getproperty(params.p, J).POC*POC)

Chlᴾ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_Chl_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :P, params)
Chlᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_Chl_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :D, params)

@inline ρᶜʰˡ(P, Chlᴾ, PARᴾ, L_day, zₘₓₗ, Feᴾ, NO₃, NH₄, PO₄, Fe, T, I, params) = 144*μ̅(P, Chlᴾ, PARᴾ, L_day, zₘₓₗ, Feᴾ, NO₃, NH₄, PO₄, Fe, I, params)/(getproperty(params.α, I)*Chlᴾ*PARᴾ/L_day)
@inline μ̅(P, Chlᴾ, PARᴾ, L_day, zₘₓₗ, Feᴾ, NO₃, NH₄, PO₄, Fe, T, I, params) = μₚ(T, I, params)*f₂(zₘₓₗ, zₑᵤ, I, params)*(1-exp(-getproperty(params.α, I)*(Chlᴾ/P)*PARᴾ/(L_day*μₚ(T, I, params)*Lₗᵢₘ(P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, I, params))))*Lₗᵢₘ(P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, I, params)
@inline Phyto_Chl_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, I, params) = (
    (1-getproperty(params.δ, I))*(120*params.θᶜʰˡₘᵢₙ+(getproperty(params.θᶜʰˡₘₐₓₓ, I) - params.θᶜʰˡₘᵢₙ)*ρᶜʰˡ((I==:P ? [P, Chlᴾ, PARᴾ] : [D, Chlᴰ, PARᴰ])..., params.L_day(t), zₘₓₗ(x, y, params), Feᴾ, NO₃, NH₄, PO₄, Fe, T, I, params))*μᴾ(params.L_day(t), zₘₓₗ(x, y, params), zₑᵤ(x, y, params), (I==:P ? [P, Chlᴾ, PARᴾ] : [D, Chlᴰ, PARᴰ])..., NO₃, NH₄, PO₄, Fe, PARᴾ, T, I, params)*(I==:P ? P : D)
    -getproperty(params.m, I)*(I==:P ? P : D)/(params.Kₘ+(I==:P ? P : D))*(I==:P ? Chlᴾ : Chlᴰ)
    -(zₘₓₗ(x, y, params)<z ? params.shₘₓₗ : params.shₛᵤ)*(I==:P ? params.wᴾ : params.wᵖ + params.wₘₐₓᴰ*(1-Lₗᵢₘ(D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, :D, params)))*(I==:P ? P*Chlᴾ : D*Chlᴰ)
    -(I==:P ? Chlᴾ/P : Chlᴰ/D)*(g(P, D, POC, I, :Z, params)*Z+g(P, D, POC, I, :M, params)*M))

Feᴾ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_Fe_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :P, params)
Feᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = Phyto_Fe_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, :D, params)

@inline μᶠᵉ(L_day, zₘₓₗ, zₑᵤ, P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, PARᴾ, T, I, params) = getproperty(params.θᶠᵉₘₐₓ, I)*Lᶠᵉ₁(Fe, Feᴾ, I)*Lᶠᵉ₂(Fe, Feᴾ, I)*μᴾ(L_day, zₘₓₗ, zₑᵤ, P, Chlᴾ, Feᴾ, NO₃, NH₄, PO₄, Fe, PARᴾ, T, I, params)
@inline θᶠᵉₘᵢₙ(P, Chlᴾ, NO₃, NH₄, I, params) = 0.0016*Chlᴾ/(55.85*P) + 1.21e-5*14*Lₙ(NO₃, NH₄, I, params)/(55.85*7.625)*1.5+1.15*14*L_NO₃(NO₃, NH₄, I, params)/(55.85*7.625) #Lₙ could be meant to be L_NH₄?
@inline Phyto_Fe_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, I, params) = (
    (1-getproperty(params.δ, I))*μᶠᵉ(params.L_day(t), zₘₓₗ(x, y, params), zₑᵤ(x, y, params), (I==:P ? [P, Chlᴾ, Feᴾ] : [D, Chlᴰ, Feᴰ])..., NO₃, NH₄, PO₄, Fe, I==:P ? PARᴾ : PARᴰ, T, I, params)*(I==:P ? P : D)
    -getproperty(params.m, I)*(I==:P ? P : D)/(params.Kₘ+(I==:P ? P : D))*(I==:P ? Feᴾ : Feᴰ)
    -(zₘₓₗ(x, y, params)<z ? params.shₘₓₗ : params.shₛᵤ)*(I==:P ? params.wᴾ : params.wᵖ + params.wₘₐₓᴰ*(1-Lₗᵢₘ(D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, :D, params)))*(I==:P ? P*Feᴾ : D*Feᴰ)
    -(I==:P ? Feᴾ/P : Feᴰ/D)*(g(P, D, POC, I, :Z, params)*Z+g(P, D, POC, I, :M, params)*M))

@inline θˢⁱᴰₒₚₜ() = params.θˢⁱᴰₘ*Lₗᵢₘˢⁱᴰ₁(Si, params)*min(5.4, (4.4*exp(-4.23*Fₗᵢₘˢⁱᴰ₁())*Fₗᵢₘˢⁱᴰ₂(Si, params)+1)*(1+2*Lₗᵢₘˢⁱᴰ₂(Si, ϕ(y, params), params)))

@inline Lₗᵢₘˢⁱᴰ₁(Si, params) = Si/(Si + params.Kₛᵢ¹)
@inline Lₗᵢₘˢⁱᴰ₂(Si, ϕ, params) = ϕ < 0 ? Si^3/(Si^3+(params.Kₛᵢ²)^3) : 0

@inline Fₗᵢₘˢⁱᴰ₁() = min(
    μᴾ(params.L_day(t), zₘₓₗ(x, y, params), zₑᵤ(x, y, params), D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, PARᴰ, T, :D, params)/(μₚ(T, :D, params)*Lₗᵢₘ(D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, :D, params)),
    Lₚₒ₄(PO₄, D, :D, params),
    Lₙ(NO₃, NH₄, :D, params),
    L_Fe(D, Chlᴰ, Feᴰ, NO₃, NH₄, :D, params)
)
@inline Fₗᵢₘˢⁱᴰ₂(Si, params) = min(1, 2.2*max(0, Lₗᵢₘˢⁱᴰ₁(Si, params)-0.5))

Siᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, NUM, Feᴾᴼ, Siᴾᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ::AbstractFloat, PARᵟ::AbstractFloat, T::AbstractFloat, params) = (
    θˢⁱᴰₒₚₜ()*(1-params.δ.D)*μᴾ(params.L_day(t), zₘₓₗ(x, y, params), zₑᵤ(x, y, params), D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, PARᴰ, T, :D, params)*D
    -params.m.D*D/(params.Kₘ+D)*Siᴰ
    -(zₘₓₗ(x, y, params)<z ? params.shₘₓₗ : params.shₛᵤ)*(params.wᵖ + params.wₘₐₓᴰ*(1-Lₗᵢₘ(D, Chlᴰ, Feᴰ, NO₃, NH₄, PO₄, Fe, :D, params)))*D*Siᴰ
    -(Siᴰ/D)*(g(P, D, POC, :D, :Z, params)*Z + g(P, D, POC, :D, :M, params)*M))


tracers=(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :NUM, :Feᴾᴼ, :Siᴾᴼ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :O₂)
optional_tracers=NamedTuple()

#forcing_functions=(P=P_forcing, D=D_forcing, Chlᴾ=Chlᴾ_forcing, Chlᴰ=Chlᴰ_forcing, Feᴾ=Feᴾ_forcing, Feᴰ=Feᴰ_forcing, Siᴰ=Siᴰ_forcing, Z=Z_forcing, M=M_forcing, DOC=DOC_forcing, POC=POC_forcing, NUM=NUM_forcing, Feᴾᴼ=Feᴾᴼ_forcing, Siᴾᴼ=Siᴾᴼ_forcing, NO₃=NO₃_forcing, NH₄=NH₄_forcing, PO₄=PO₄_forcing, Fe=Fe_forcing, Si=Si_forcing, CaCO₃=CaCO₃_forcing, DIC=DIC_forcing, O₂=O₂_forcing)
#sinking=(POC=POC_sinking, NUM=NUM_sinking)

#parameters
const defaults = (
    
)
#end # module