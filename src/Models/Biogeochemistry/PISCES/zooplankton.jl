#grazing of zooplankton, here food_consumption = Pᵢᶻ*(I - Iₜₕᵣᶻ), food_total = sum(pⱼᶻ*J) and other variables are as eq 26a
#foods are (; P, D, POC, Z) for mezozooplankton and  (; P, D, POC) for microzooplankton
@inline g(food_consumption, gₘ, Fₗᵢₘ, F, K_G, food_total) = gₘ*(Fₗᵢₘ/(F+eps(0.0)))*max(0, food_consumption/(K_G+food_total))
@inline Fᵢ(p, Jₜₕᵣ, C) = p*max(0, C - Jₜₕᵣ)
@inline F(p, Jₜₕᵣ, prey) = sum(Fᵢ(p[I], Jₜₕᵣ, C) for (I, C) in pairs(prey))
@inline Fₗᵢₘ(p, Jₜₕᵣ, Fₜₕᵣ, prey) = max(0, F(p, Jₜₕᵣ, prey)-min(0.5*F(p, Jₜₕᵣ, prey), Fₜₕᵣ))
@inline food_total(p, prey) = sum([@inbounds p[I]*C for (I, C) in pairs(prey)])

@inline function θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴾ, T, zₘₓₗ, zₑᵤ, K_PO₄, K_NO₃, K_NH₄, θᶠᵉₒₚₜ, Kₛᵢᴰᵐⁱⁿ, Si̅, Kₛᵢ, μₘₐₓ⁰, b, t_dark, L_day, α)
    #Not convinced this is correct as no explicit definition of θᴺᴵ after eq 27b so infered from NEMO source code
    Lₚₒ₄ᴾ = L_mondo(PO₄, K_PO₄)
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, K_NO₃, K_NH₄)
    Lₙᴾ = Lₙₒ₃ᴾ + L_NH₄(NO₃, NH₄, K_NO₃, K_NH₄) #eq 6c
    L_Feᴾ = L_Fe(P, Chlᴾ, Feᴾ, θᶠᵉₒₚₜ, Lₙᴾ, Lₙₒ₃ᴾ)
    Lₛᵢᴾ = L_mondo(Si, Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2/(Kₛᵢ^2+Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴾ = min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ, Lₛᵢᴾ)

    μₚ = μₘₐₓ⁰*b^T

    return f₂(zₘₓₗ, zₑᵤ, t_dark)*Cₚᵣₒ(P, Chlᴾ, PARᴾ, L_day, α, μₚ, Lₗᵢₘᴾ)*min(Lₚₒ₄ᴾ, Lₙᴾ)
end

function Z_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #grazing
    prey =  (; P, D, POC)

    L_day = params.L_day(t)
    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    Σgᶻ = gᶻₚ + gᶻ_D + gᶻₚₒ

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #total Fe and N grazing
    gᶻ_Fe = gᶻₚ*Feᴾ/(P + eps(0.0))+ gᶻ_D*Feᴰ/(D + eps(0.0))+ gᶻₚₒ*Feᴾᴼ/(POC + eps(0.0)) 
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    #efficiency (from food quality)
    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*Σgᶻ + eps(0.0)), gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))
    eᶻ = eₙᶻ*min(params.eₘₐₓ.Z, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))

    #temperature dependency
    f_z = params.b.Z^T

    #getting grazed
    gᴹ_z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.M), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, food_total(params.p.M, (; P, D, POC, Z)))

    return eᶻ*Σgᶻ*Z - gᴹ_z*Z - params.m.Z*Z^2*f_z - params.r.Z*f_z*(Z/(params.Kₘ+Z) + 3*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))*Z
end

function M_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #grazing
    prey =  (; P, D, POC, Z)

    L_day = params.L_day(t)
    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_Z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ)
    Σgᴹ = gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹ_Z

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #temperature dependency
    f_M = params.b.M^T

    #flux feeding on POC - will need to be changed for Kriest Model
    gᴹₚₒ_FF = params.g_FF*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_FF*f_M*(params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    #total Fe and N grazing - based on paragraph after eq 32 I don't think flux feeding is included in the quality
    gᴹ_Fe = gᴹₚ*Feᴾ/(P+eps(0.0))+ gᴹ_D*Feᴰ/(D+eps(0.0))+ gᴹₚₒ*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z)*params.θᴺᶜ

    #efficiency (from food quality)
    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ+eps(0.0)), gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ+eps(0.0)))
    eᴹ = eₙᴹ*min(params.eₘₐₓ.M, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ+eps(0.0)))

    #Mortality from upper trophic feeding or density dependency (e.g. viral disease) - if a separate feeding model is coupled fallback to microzooplankton mortality
    #Also need to modify POC since if upper trophic feeding on need to add this mortality to POC (like microzooplankton) 
    mᴹ = ifelse(params.upper_trophic_feeding, params.m.M, params.m.Z)

    return eᴹ*(Σgᴹ + gᴹₚₒ_FF + gᴹ_GO_FF)*M - mᴹ*M^2*f_M - params.r.M*f_M*(M/(params.Kₘ+M) + 3*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))*M
end