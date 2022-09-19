#grazing of zooplankton, here food_consumption = Pᵢᶻ*(I - Iₜₕᵣᶻ), food_total = sum(pⱼᶻ*J) and other variables are as eq 26a
#foods are (; P, D, POC, Z) for mezozooplankton and  (; P, D, POC) for microzooplankton
@inline g(food_consumption, gₘ, Fₗᵢₘ, F, K_G, food_total) = gₘ*(Fₗᵢₘ/F)*max(0, food_consumption/(K_G+food_total)
@inline F(p, Jₜₕᵣ, prey) = sum([p[I]*max(0, C-Jₜₕᵣ[I]) for (I, C) in pairs(prey)])
@inline Fₗᵢₘ(p, Jₜₕᵣ, Fₜₕᵣ, prey) = max(0, F(p, Jₜₕᵣ, prey)-min(0.5*F(p, Jₜₕᵣ, prey), Fₜₕᵣ))
@inline food_total(p, prey) = sum([@inbounds p[I]*C for (I, C) in pairs(prey)])

function Z_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾᴼ, Siᴳᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ, PARᴰ, T, zₘₓₗ, zₑᵤ, ϕ, params)
    #grazing
    prey =  (; P, D, POC)

    L_day = params.L_day(t)
    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(params.p.Z.P*(P-params.Jₜₕᵣ.Z.P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(params.p.Z.D*(D-params.Jₜₕᵣ.Z.D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z.POC), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    Lₚₒ₄ᴾ = L_mondo(PO₄, K(P, params.Kᵐⁱⁿ.PO₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P))
    L_NO₃ᴾ = L_NO₃ᴾ(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P) #eq 6c
    L_Feᴾ = L_Fe(P, Chlᴾ, Feᴾ, params.θᶠᵉₒₚₜ.P, Lₙᴾ, Lₙₒ₃ᴾ)
    Lₗᵢₘᴾ = min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ)
    #Not convinced this is correct as no explicit definition of θᴺᴵ after eq 27b so infered from NEMO source code
    θᴺᴾ = f₂(zₘₓₗ, zₑᵤ, params.t_dark.P)*Cₚᵣₒ(P, Chlᴾ, PARᴾ, L_day, params.α.P, μₚ, Lₗᵢₘᴾ)*min(Lₚₒ₄ᴾ, Lₙᵖ)

    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, params.Kᵐⁱⁿ.PO₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D))
    L_NO₃ᴰ = L_NO₃ᴾ(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D)
    Lₙᴰ = L_NO₃ᴰ + L_NH₄(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D) #eq 6c
    L_Feᴰ = L_Fe(P, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)
    Lₛᵢᴰ = L_mondo(Si, params.Kₛᵢᴰᵐⁱⁿ + 7*params.Si̅^2/(params.Kₛᵢ^2+params.Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴰ = min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ)
    θᴺᴰ = f₂(zₘₓₗ, zₑᵤ, params.t_dark.D)*Cₚᵣₒ(D, Chlᴰ, PARᴰ, L_day, params.α.D, μ_d, Lₗᵢₘᴰ)*min(Lₚₒ₄ᴰ, Lₙᴰ)

    #total Fe and N grazing
    gᶻ_Fe = gᶻₚ*Feᴾ/P + gᶻ_D*Feᴰ/D + gᶻₚₒ*Feᴾᴼ/POC
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    #efficiency (from food quality)
    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*(gᶻₚ + gᶻ_D + gᶻₚₒ)), gᶻ_Fe/(params.θᶠᵉ.Z*(gᶻₚ + gᶻ_D + gᶻₚₒ)))
    eᶻ = eₙᶻ*min(params.eₘₐₓᶻ, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*(gᶻₚ + gᶻ_D + gᶻₚₒ)))

    #temperature dependency
    f_z = params.b.Z^T

    #getting grazed
    gᴹ_z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.Z), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, food_total(params.p.M, (; P, D, POC, Z)))

    return eᶻ*(gᶻₚ + gᶻ_D + gᶻₚₒ)*Z - gᴹ_z *M - params.m.Z*Z^2*f_z - params.r.Z*f_z*(Z/(params.Kₘ+z) + 3*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))*Z
end


function M_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾᴼ, Siᴳᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ, PARᴰ, T, zₘₓₗ, zₑᵤ, ϕ, params)
    #grazing
    prey =  (; P, D, POC, Z)

    L_day = params.L_day(t)
    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚ = g(params.p.M.P*(P-params.Jₜₕᵣ.M.P), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_D = g(params.p.M.D*(D-params.Jₜₕᵣ.M.D), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M.POC), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 

    Lₚₒ₄ᴾ = L_mondo(PO₄, K(P, params.Kᵐⁱⁿ.PO₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P))
    L_NO₃ᴾ = L_NO₃ᴾ(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P) #eq 6c
    L_Feᴾ = L_Fe(P, Chlᴾ, Feᴾ, params.θᶠᵉₒₚₜ.P, Lₙᴾ, Lₙₒ₃ᴾ)
    Lₗᵢₘᴾ = min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ)
    #Not convinced this is correct as no explicit definition of θᴺᴵ after eq 27b so infered from NEMO source code
    θᴺᴾ = f₂(zₘₓₗ, zₑᵤ, params.t_dark.P)*Cₚᵣₒ(P, Chlᴾ, PARᴾ, L_day, params.α.P, μₚ, Lₗᵢₘᴾ)*min(Lₚₒ₄ᴾ, Lₙᵖ)

    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, params.Kᵐⁱⁿ.PO₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D))
    L_NO₃ᴰ = L_NO₃ᴾ(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D)
    Lₙᴰ = L_NO₃ᴰ + L_NH₄(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D) #eq 6c
    L_Feᴰ = L_Fe(P, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)
    Lₛᵢᴰ = L_mondo(Si, params.Kₛᵢᴰᵐⁱⁿ + 7*params.Si̅^2/(params.Kₛᵢ^2+params.Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴰ = min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ)
    θᴺᴰ = f₂(zₘₓₗ, zₑᵤ, params.t_dark.D)*Cₚᵣₒ(D, Chlᴰ, PARᴰ, L_day, params.α.D, μ_d, Lₗᵢₘᴰ)*min(Lₚₒ₄ᴰ, Lₙᴰ)

    #temperature dependency
    f_M = params.b.M^T

    #flux feeding on POC - will need to be changed for Kriest Model
    gᴹₚₒ_FF = params.g_ff*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_ff*f_M*(params._GOᵐⁱⁿ + (200-params._GOᵐⁱⁿ)*max(0, z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    #total Fe and N grazing
    gᴹ_Fe = gᴹₚ*Feᴾ/P + gᴹ_D*Feᴰ/D + (gᴹₚₒ+gᴹₚₒ_FF)*Feᴾᴼ/POC + gᴹ_GO_FF*Feᴳᴼ/GOC
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ+gᴹₚₒ_FF+gᴹ_GO_FF)*params.θᴺᶜ 

    #efficiency (from food quality)
    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*(gᴹₚ + gᴹ_D + gᴹₚₒ+gᴹₚₒ_FF+gᴹ_GO_FF), gᴹ_Fe/(params.θᶠᵉ.M*(gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹₚₒ_FF + gᴹ_GO_FF)))
    eᴹ = eₙᴹ*min(params.eₘₐₓᴹ, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*(gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹₚₒ_FF + gᴹ_GO_FF)))

    #Mortality from upper trophic feeding or density dependency (e.g. viral disease) - if a separate feeding model is coupled fallback to microzooplankton mortality
    #Also need to modify POC since if upper trophic feeding on need to add this mortality to POC (like microzooplankton) 
    mᴹ = ifelse(params.upper_trophic_feeding, params.mᴹ, params.mᶻ)

    return eᴹ*(gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹₚₒ_FF + gᴹ_GO_FF)*M - mᴹ*M^2*f_M - params.r.M*f_M*(M/(params.Kₘ+M) + 3*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))*M
end