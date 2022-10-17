function DIC_forcing(i, j, k, grid, clock, model_fields, params) #needs to be discrete because of \bar{PAR}
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)
    L_day = params.L_day(t)

    #phytoplankton growth
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #microzooplankton grazing
    prey =  (; P, D, POC)
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

    gᶻ_Fe = gᶻₚ*Feᴾ/(P+eps(0.0))+ gᶻ_D*Feᴰ/(D+eps(0.0))+ gᶻₚₒ*Feᴾᴼ/(POC+eps(0.0)) 
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*Σgᶻ + eps(0.0)), gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))
    eᶻ = eₙᶻ*min(params.eₘₐₓ.Z, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))

    #Macrozooplankton grazing
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
    
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    Lₙᴾ = Lₙₒ₃ᴾ + L_NH₄(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)) 
    r_CaCO₃ = R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    f_M = params.b.M^T

    gᴹₚₒ_FF = params.g_FF*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_FF*f_M*(params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    gᴹ_Fe = gᴹₚ*Feᴾ/(P+eps(0.0))+ gᴹ_D*Feᴰ/(D+eps(0.0))+ gᴹₚₒ*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z + eps(0.0))*params.θᴺᶜ

    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ + eps(0.0)), gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ + eps(0.0)))
    eᴹ = eₙᴹ*min(params.eₘₐₓ.M, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ + eps(0.0)))

    #Upper trophic feeding waste
    if params.upper_trophic_feeding
        Rᵤₚᴹ = (1 - params.σ.M - params.eₘₐₓ.M)*(1/(1-params.eₘₐₓ.M))*params.m.M*f_M*M^2
    else
        Rᵤₚᴹ = 0.0
    end

    #Bacterial action
    Bactᵢⱼₖ = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))
    ΔO₂ᵢⱼₖ = ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)

    bFe = Fe
    Lᴮᵃᶜᵗ = Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.Kₙₒ₃.Bact, params.Kₙₕ₄.Bact)*L_mondo(DOC, params.K_DOC)
    μ_Bact = params.λ_DOC*params.bₚ*Bactᵢⱼₖ*DOC*Lᴮᵃᶜᵗ#eq 33a, b say Lₗᵢₘᴮᵃᶜᵗ but think this must be Lᴮᵃᶜᵗ

    Remin = min(O₂/params.O₂ᵘᵗ, μ_Bact*(1-ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)))
    Denit = min(NO₃/params.rₙₒ₃, μ_Bact*ΔO₂ᵢⱼₖ)
    
    #dissolusion - this is not fully described in PISCES paper so copied from NEMO's code and Mucci 1983 10.2475/ajs.283.7.780
    Conc_Ca = 1.0287e-2 * S / 35.16504 #zcalcon
    logKɔ = -171.945 - 0.077993*T +2903.293/T + 71.595*log(10, T)
    B = -0.77712 + 2.8426*T + 178.34/T 
    Ksp = 10^(logKɔ + B*S^0.5 + -0.07711*S +4.1249*S^1.5)
    Ω_Ca = Conc_Ca*CaCO₃/Ksp
    #CO₃⁻²ₛₐₜ = Ksp / Conc_Ca

    ΔCO₃⁻² = max(0, 1 - Ω_Ca) 

    λ_CaCO₃ = ifelse(Ω_Ca < 0.8, params.λ_CaCO₃*ΔCO₃⁻²^params.nca, params.λ_CaCO₃*(0.2^(params.nca - 0.11))*ΔCO₃⁻²^0.11)  
    
    #Production
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 
    
    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)
    
    p_waste = params.m.P*P*L_mondo(P, params.Kₘ)+sh*params.wᴾ*P^2

    return (
        params.γ.Z*(1 - eᶻ - params.σ.Z)*Σgᶻ*Z 
        + params.γ.M*(1 - eᴹ - params.σ.M)*(Σgᴹ + gᴹₚₒ_FF + gᴹ_GO_FF)*M 
        + params.γ.M*Rᵤₚᴹ + Remin + Denit
        - r_CaCO₃*(params.η.Z*gᶻₚ*Z + params.η.M*gᴹₚ*M + 0.5*p_waste)
        + λ_CaCO₃*CaCO₃ 
        - μᴰ*D - μᴾ*P 
    ) 
end

function Alk_forcing(i, j, k, grid, clock, model_fields, params) #needs to be discrete because of \bar{PAR}
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)
    L_day = params.L_day(t)

    #phytoplankton growth
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    L_NH₄ᴾ =  + L_NH₄(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    μₙₒ₃ᴾ = μᴾ * L_NO₃ᴾ / (L_NO₃ᴾ + L_NH₄ᴾ + eps(0.0))
    μₙₕ₄ᴾ = μᴾ * L_NH₄ᴾ / (L_NO₃ᴾ + L_NH₄ᴾ + eps(0.0))

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    L_NO₃ᴰ = L_NO₃(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    L_NH₄ᴰ =  + L_NH₄(NO₃, NH₄, K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    μₙₒ₃ᴰ = μᴰ * L_NO₃ᴰ / (L_NO₃ᴰ + L_NH₄ᴰ + eps(0.0))
    μₙₕ₄ᴰ = μᴰ * L_NH₄ᴰ / (L_NO₃ᴰ + L_NH₄ᴰ + eps(0.0))

    #microzooplankton grazing
    prey =  (; P, D, POC)
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
    
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    Lₙᴾ = Lₙₒ₃ᴾ + L_NH₄(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)) 
    r_CaCO₃ = R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    gᶻ_Fe = gᶻₚ*Feᴾ/(P+eps(0.0))+ gᶻ_D*Feᴰ/(D+eps(0.0))+ gᶻₚₒ*Feᴾᴼ/(POC+eps(0.0)) 
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*Σgᶻ + eps(0.0)), gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))
    eᶻ = eₙᶻ*min(params.eₘₐₓ.Z, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))

    #Macrozooplankton grazing
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

    f_M = params.b.M^T

    gᴹₚₒ_FF = params.g_FF*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_FF*f_M*(params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    gᴹ_Fe = gᴹₚ*Feᴾ/(P+eps(0.0))+ gᴹ_D*Feᴰ/(D+eps(0.0))+ gᴹₚₒ*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z)*params.θᴺᶜ

    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ + eps(0.0)), gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ + eps(0.0)))
    eᴹ = eₙᴹ*min(params.eₘₐₓ.M, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ + eps(0.0)))

    #Upper trophic feeding waste
    if params.upper_trophic_feeding
        Rᵤₚᴹ = (1 - params.σ.M - params.eₘₐₓ.M)*(1/(1-params.eₘₐₓ.M))*params.m.M*f_M*M^2
    else
        Rᵤₚᴹ = 0.0
    end

    #Bacterial action
    Bactᵢⱼₖ = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))
    ΔO₂ᵢⱼₖ = ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)

    bFe = Fe
    Lᴮᵃᶜᵗ = Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.Kₙₒ₃.Bact, params.Kₙₕ₄.Bact)*L_mondo(DOC, params.K_DOC)
    μ_Bact = params.λ_DOC*params.bₚ*Bactᵢⱼₖ*DOC*Lᴮᵃᶜᵗ#eq 33a, b say Lₗᵢₘᴮᵃᶜᵗ but think this must be Lᴮᵃᶜᵗ

    Remin = min(O₂/params.O₂ᵘᵗ, μ_Bact*(1-ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)))
    Denit = min(NO₃/params.rₙₒ₃, μ_Bact*ΔO₂ᵢⱼₖ)

    kₘₓₗ = round(Int, fractional_z_index(-zₘₓₗ, Center(), grid))+1
    PAR̄ = (sum(model_fields.PAR¹[i, j, kₘₓₗ:grid.Nz].*grid.zᵃᵃᶜ[kₘₓₗ:grid.Nz]) + sum(model_fields.PAR²[i, j, kₘₓₗ:grid.Nz].*grid.zᵃᵃᶜ[kₘₓₗ:grid.Nz]) + sum(model_fields.PAR³[i, j, kₘₓₗ:grid.Nz].*grid.zᵃᵃᶜ[kₘₓₗ:grid.Nz]))/zₘₓₗ
    Nitrif =  params.λₙₕ₄*NH₄*(1-ΔO₂ᵢⱼₖ)/(1+PAR̄)

    μₚ = params.μₘₐₓ⁰*params.bₚ^T
    bFe = Fe

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P

    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ))
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)) 
    Lₙᴰᶻ = ifelse(Lₙᴾ >= 0.8, 0.01, 1 - Lₙᴾ)
    N_fix = params.N_fixᵐ*max(0, μₚ - 2.15)*Lₙᴰᶻ*min(L_mondo(bFe, params.K_Feᴰᶻ), L_mondo(PO₄, Kₚₒ₄ᵐⁱⁿ))*(1 - exp(-PAR/params.E_fix))
    
    #dissolusion - this is not fully described in PISCES paper so copied from NEMO's code and Mucci 1983 10.2475/ajs.283.7.780
    Conc_Ca = 1.0287e-2 * S / 35.16504 #zcalcon
    logKɔ = -171.945 - 0.077993*T +2903.293/T + 71.595*log(10, T)
    B = -0.77712 + 2.8426*T + 178.34/T 
    Ksp = 10^(logKɔ + B*S^0.5 + -0.07711*S +4.1249*S^1.5)
    Ω_Ca = Conc_Ca*CaCO₃/Ksp
    #CO₃⁻²ₛₐₜ = Ksp / Conc_Ca

    ΔCO₃⁻² = max(0, 1 - Ω_Ca)

    λ_CaCO₃ = ifelse(Ω_Ca < 0.8, params.λ_CaCO₃*ΔCO₃⁻²^params.nca, params.λ_CaCO₃*(0.2^(params.nca - 0.11))*ΔCO₃⁻²^0.11) 
    
    #Production
    z_food_total = food_total(params.p.Z, (; P, D, POC)) #could store this as an auxiliary field that is forced so it doesn't need to be computed for P, D, POC and Z
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) #may be inconcenient/impossible to get it to calculate it first though, could make this descrete 
    
    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)
    
    p_waste = params.m.P*P*L_mondo(P, params.Kₘ)+sh*params.wᴾ*P^2

    return (
        params.θᴺᶜ*params.γ.Z*(1 - eᶻ - params.σ.Z)*Σgᶻ*Z 
        + params.θᴺᶜ*params.γ.M*(1 - eᴹ - params.σ.M)*(Σgᴹ + gᴹₚₒ_FF + gᴹ_GO_FF)*M #must be a typo on line 4 of eq 81 since Rᵤₚᴹ isn't a per M value and shouldn't get multipled by θ and γ twice? 
        + params.θᴺᶜ*params.γ.M*Rᵤₚᴹ 
        + params.θᴺᶜ*(Remin + (1 + params.rₙₒ₃)*Denit)
        - 2*r_CaCO₃*(params.η.Z*gᶻₚ*Z + params.η.M*gᴹₚ*M + 0.5*p_waste)
        + 2*λ_CaCO₃*CaCO₃ 
        + params.θᴺᶜ*(N_fix - 2*Nitrif)
        + params.θᴺᶜ*ΔO₂ᵢⱼₖ*(params.rₙₕ₄ - 1)*params.λₙₕ₄*NH₄
        + params.θᴺᶜ*((μₙₒ₃ᴾ - μₙₕ₄ᴾ)*P + (μₙₒ₃ᴰ - μₙₕ₄ᴰ)*D)
    ) 
end