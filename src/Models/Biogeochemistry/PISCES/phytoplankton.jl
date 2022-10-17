@inline Cₚᵣₒ(P, Chlᴾ, PARᴾ, L_day, α, μₚ, Lₗᵢₘᴾ)=1-exp(-α*(Chlᴾ/(P+eps(0.0)))*PARᴾ/(L_day*μₚ*Lₗᵢₘᴾ + eps(0.0)))

@inline f₁(L_day) = 1.5*L_day/(0.5+L_day)#eq 3a
@inline t_dark(zₘₓₗ, zₑᵤ) = max(0, zₘₓₗ-zₑᵤ)^2/86400#eq 3b,c
@inline f₂(zₘₓₗ, zₑᵤ, t_dark_lim) = 1 - t_dark(zₘₓₗ, zₑᵤ)/(t_dark(zₘₓₗ, zₑᵤ)+t_dark_lim) #eq 3d

@inline L_NO₃(NO₃, NH₄, K_NO₃, K_NH₄) = K_NO₃*NH₄/(K_NO₃*K_NH₄+K_NH₄*NO₃+K_NO₃*NH₄) #eq 6e
@inline L_NH₄(NO₃, NH₄, K_NO₃, K_NH₄) = K_NH₄*NO₃/(K_NO₃*K_NH₄+K_NH₄*NO₃+K_NO₃*NH₄) #eq 6d
#@inline Lₙ(NO₃, NH₄, K_NO₃, K_NH₄) = L_NO₃(NO₃, NH₄, K_NO₃, K_NH₄) + L_NH₄(NO₃, NH₄, K_NO₃, K_NH₄)
@inline L_mondo(C, K) = C/(eps(0.0) + C + K) #eq 6b

@inline L_Fe(P, Chlᴾ, Feᴾ, θᶠᵉₒₚₜ, Lₙᴾ, Lₙₒ₃ᴾ) = min(1, max(0, ((Feᴾ/(P+eps(0.0))) - θᶠᵉₘᵢₙ(P, Chlᴾ, Lₙᴾ, Lₙₒ₃ᴾ))/θᶠᵉₒₚₜ)) #eq 6f
@inline θᶠᵉₘᵢₙ(P, Chlᴾ, Lₙᴾ, Lₙₒ₃ᴾ) = 0.0016*Chlᴾ/(55.85*(P+eps(0.0))) + 1.21e-5*14*Lₙᴾ/(55.85*7.625)*1.5+1.15*14*Lₙₒ₃ᴾ/(55.85*7.625) #eq 20 -> Lₙ could be meant to be L_NH₄?

@inline K(P, Kᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ) = Kᵐⁱⁿ*(P₁(P, Pₘₐₓ)+Sᵣₐₜ*P₂(P, Pₘₐₓ))/(P₁(P, Pₘₐₓ)+P₂(P, Pₘₐₓ))
@inline P₁(P, Pₘₐₓ) = min(P, Pₘₐₓ)
@inline P₂(P, Pₘₐₓ) = min(0, P - Pₘₐₓ)

@inline function μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴾ, T, zₘₓₗ, zₑᵤ, K_PO₄, K_NO₃, K_NH₄, θᶠᵉₒₚₜ, Kₛᵢᴰᵐⁱⁿ, Si̅, Kₛᵢ, μₘₐₓ⁰, b, t_dark, L_day, α)
    #eq2a and 11b
    #Not convinced this is correct as no explicit definition of θᴺᴵ after eq 27b so infered from NEMO source code
    Lₚₒ₄ᴾ = L_mondo(PO₄, K_PO₄)
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, K_NO₃, K_NH₄)
    Lₙᴾ = Lₙₒ₃ᴾ + L_NH₄(NO₃, NH₄, K_NO₃, K_NH₄) #eq 6c
    L_Feᴾ = L_Fe(P, Chlᴾ, Feᴾ, θᶠᵉₒₚₜ, Lₙᴾ, Lₙₒ₃ᴾ)
    Lₛᵢᴾ = L_mondo(Si, Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2/(Kₛᵢ^2+Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴾ = min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ, Lₛᵢᴾ)

    μₚ = μₘₐₓ⁰*b^T

    return  μₚ*f₁(L_day)*f₂(zₘₓₗ, zₑᵤ, t_dark)*Cₚᵣₒ(P, Chlᴾ, PARᴾ, L_day, α, μₚ, Lₗᵢₘᴾ)*min(Lₚₒ₄ᴾ, Lₙᴾ), Lₗᵢₘᴾ
end

#perhaps should add an auxiliary forcing before to compute all of the reused values such as z_food_total, F_lim^Z etc.? which otherwise get computed 4 times minimum
#think the trade off here varies for CPU vs GPU where we might want to not recompute them on CPU but we may want to store less in memory on GPU
function P_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)

    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.P)*μᴾ*P - params.m.P*P^2/(params.Kₘ+P) - sh*params.wᴾ*P^2 - gᶻₚ*Z - gᴹₚ*M#eq 1
end

function D_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.D)*μᴰ*D - params.m.D*D^2/(params.Kₘ+D) - sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D^2 - gᶻ_D*Z - gᴹ_D*M#eq 9
end

function Chlᴾ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)

    μ̌ᴾ = (μᴾ/f₁(L_day))
    ρᶜʰˡᴾ = 144*μ̌ᴾ*P/(eps(0.0) + params.α.P*Chlᴾ*PARᴾ/L_day)
    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.P)*(12*params.θᶜʰˡₘᵢₙ + (params.θᶜʰˡₘₐₓ.P-params.θᶜʰˡₘᵢₙ)*ρᶜʰˡᴾ)*μᴾ*P - params.m.P*P*Chlᴾ/(params.Kₘ+P) - sh*params.wᴾ*P*Chlᴾ - (Chlᴾ/(P+eps(0.0)))*(gᶻₚ*Z + gᴹₚ*M)#eq 14
end

function Chlᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    μ̌ᴰ = (μᴰ/f₁(L_day)) #eq 15a
    ρᶜʰˡᴰ = 144*μ̌ᴰ*D/(eps(0.0) + params.α.D*Chlᴰ*PARᴰ/L_day) #eq 15a
    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.D)*(12*params.θᶜʰˡₘᵢₙ + (params.θᶜʰˡₘₐₓ.D-params.θᶜʰˡₘᵢₙ)*ρᶜʰˡᴰ)*μᴰ*D - params.m.D*D*Chlᴰ/(params.Kₘ+D) - sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D*Chlᴰ - (Chlᴰ/(D+eps(0.0)))*(gᶻ_D*Z + gᴹ_D*M)#eq 14
end

function Feᴾ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    μₚ = (params.μₘₐₓ⁰*params.bₚ^T)

    P₂ = max(0, P - params.Pₘₐₓ.P)
    P₁ = P - P₂
    K_Feᶠᵉᴾ = params.K_Feᵐⁱⁿ.P*(P₁ + params.Sᵣₐₜ.P*P₂)/(eps(0.0) + P₁ + P₂)

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P

    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    Lₙᴾ = Lₙₒ₃ᴾ + L_NH₄(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)) #eq 6c
    L_Feᴾ = L_Fe(P, Chlᴾ, Feᴾ, params.θᶠᵉₒₚₜ.P, Lₙᴾ, Lₙₒ₃ᴾ)

    bFe = Fe #with the complex chemsitry implimented this would not just be Fe
    μᶠᵉᴾ = params.θᶠᵉₘₐₓ.P*L_mondo(bFe, K_Feᶠᵉᴾ)*((4-4.5*L_Feᴾ)/(L_Feᴾ+0.5))*((1 - (Feᴾ/(P+eps(0.0)))/params.θᶠᵉₘₐₓ.P)/(1.05 - (Feᴾ/(P+eps(0.0)))/params.θᶠᵉₘₐₓ.P))*μₚ

    #Not really sure why it defined θᶠᵉₘᵢₙ, perhaps a typo or used somewhere else
   
    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.P)*μᶠᵉᴾ*P - params.m.P*P*Feᴾ/(params.Kₘ+P) - sh*params.wᴾ*P*Feᴾ - (Feᴾ/(P+eps(0.0)))*(gᶻₚ*Z + gᴹₚ*M)#eq 16
end

function Feᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    μ_D = (params.μₘₐₓ⁰*params.bₚ^T)

    D₂ = max(0, D - params.Pₘₐₓ.D)
    D₁ = D - D₂
    K_Feᶠᵉᴰ = params.K_Feᵐⁱⁿ.D*(D₁ + params.Sᵣₐₜ.D*D₂)/(eps(0.0) + D₁ + D₂)

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D

    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Lₙₒ₃ᴰ = L_NO₃(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Lₙᴰ = Lₙₒ₃ᴰ + L_NH₄(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ)) #eq 6c
    L_Feᴰ = L_Fe(P, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)
    Lₛᵢᴰ = L_mondo(Si, params.Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2/(params.Kₛᵢ^2+Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴰ = min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ)

    bFe = Fe #with the complex chemsitry implimented this would not just be Fe
    μᶠᵉᴰ = params.θᶠᵉₘₐₓ.D*L_mondo(bFe, K_Feᶠᵉᴰ)*((4-4.5*L_Feᴰ)/(L_Feᴰ+0.5))*((1 - (Feᴰ/(D+eps(0.0)))/params.θᶠᵉₘₐₓ.D)/(1.05 - (Feᴰ/(D+eps(0.0)))/params.θᶠᵉₘₐₓ.D))*μ_D
    
    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return (1 - params.δ.D)*μᶠᵉᴰ*D - params.m.D*D*Feᴰ/(params.Kₘ+D) - sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D*Feᴰ - (Feᴰ/(D+eps(0.0)))*(gᶻ_D*Z + gᴹ_D*M)#eq 17
end

function Siᴰ_forcing(x, y, z, t, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust, params)
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    #growth
    L_day = params.L_day(t)
    
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)
    μ_D = params.μₘₐₓ⁰*params.bₚ^T
    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Lₙₒ₃ᴰ = L_NO₃(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Lₙᴰ = Lₙₒ₃ᴰ + L_NH₄(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ)) #eq 6c
    L_Feᴰ = L_Fe(P, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)

    #optimum quota
    Lₗᵢₘ₁ˢⁱᴰ = L_mondo(Si, params.Kₛᵢ¹) #eq 23c
    Lₗᵢₘ₂ˢⁱᴰ = params.ϕ < 0 ?  L_mond(Si^3, params.Kₛᵢ²^3) : 0 #eq 23d
    Fₗᵢₘ₁ˢⁱᴰ = min(μᴰ/(eps(0.0) + μ_D*Lₗᵢₘᴰ), Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ) #eq 23a
    Fₗᵢₘ₂ˢⁱᴰ = min(1, 2.2*max(0, Lₗᵢₘ₁ˢⁱᴰ - 0.5)) #eq 23b

    θₒₚₜˢⁱᴰ = params.θˢⁱᴰₘ*Lₗᵢₘ₁ˢⁱᴰ*min(5.4, (4.4*exp(-4.23*Fₗᵢₘ₁ˢⁱᴰ)*Fₗᵢₘ₂ˢⁱᴰ+1)*(1+2*Lₗᵢₘ₂ˢⁱᴰ))#eq 22
    
    #shear
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 

    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)

    return θₒₚₜˢⁱᴰ*(1 - params.δ.D)*μᴰ*D - params.m.D*D*Siᴰ/(params.Kₘ+D) - sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D*Siᴰ - (Siᴰ/(D+eps(0.0)))*(gᶻ_D*Z + gᴹ_D*M)#eq 17
end