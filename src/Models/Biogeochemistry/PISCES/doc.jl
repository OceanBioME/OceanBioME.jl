#since Bact needs access to the whole column of Z and M need to use discrete forcings anything that depends on it
function DOC_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)

    L_day = params.L_day(t)
    #Plankton quality
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #Microzooplankton efficiency
    prey =  (; P, D, POC)

    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    Σgᶻ = gᶻₚ + gᶻ_D + gᶻₚₒ

    gᶻ_Fe = gᶻₚ*Feᴾ/(P+eps(0.0))+ gᶻ_D*Feᴰ/(D+eps(0.0))+ gᶻₚₒ*Feᴾᴼ/(POC+eps(0.0)) 
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*Σgᶻ + eps(0.0)), gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))
    eᶻ = eₙᶻ*min(params.eₘₐₓ.Z, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ + eps(0.0)))

    #Mezozooplankton efficiency
    prey =  (; P, D, POC, Z)

    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚ = g(Fᵢ(params.p.M.P, params.Jₜₕᵣ.M, P), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_Z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ)

    f_M = params.b.M^T
    
    Σgᴹ = gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹ_Z

    gᴹ_Fe = gᴹₚ*Feᴾ/(P+eps(0.0))+ gᴹ_D*Feᴰ/(D+eps(0.0))+ gᴹₚₒ*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z)*params.θᴺᶜ

    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ+ eps(0.0)), gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ+ eps(0.0)))
    eᴹ = eₙᴹ*min(params.eₘₐₓ.M, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ+ eps(0.0)))

    gᴹ_GO_FF = params.g_FF*f_M*(params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, - z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    #plankton growth
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #Upper trophic mesozooplankton feeding
    if params.upper_trophic_feeding
        Rᵤₚᴹ = (1 - params.eₘₐₓ.M- params.σ.M)*(1/(1-params.eₘₐₓ.M))*params.m.M*f_M*M^2
    else
        Rᵤₚᴹ = 0.0
    end

    #Bacterial action
    bFe = Fe
    
    Lᴮᵃᶜᵗ = Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.Kₙₒ₃.Bact, params.Kₙₕ₄.Bact)*L_mondo(DOC, params.K_DOC)
    Bactᵢⱼₖ = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))
    μ_Bact = params.λ_DOC*params.bₚ*Bactᵢⱼₖ*DOC*Lᴮᵃᶜᵗ #Bact_ref is 1 from SM of Aumount and Bopp 2006

    Remin = min(O₂/params.O₂ᵘᵗ, μ_Bact*(1-ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)))
    Denit = min(NO₃/params.rₙₒ₃, μ_Bact*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))

    #aggregation
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)
    
    ϕ₁ᴰᴼᶜ = sh*(params.a₁*DOC + params.a₂*POC)*DOC 
    ϕ₂ᴰᴼᶜ = sh*params.a₃*GOC*DOC
    ϕ₃ᴰᴼᶜ = (params.a₄*POC + params.a₅*DOC)*DOC
    return (1-params.γ.Z)*(1 - eᶻ - params.σ.Z)*Σgᶻ*Z + (1-params.γ.M)*(1 - eᴹ - params.σ.M)*(Σgᴹ + gᴹ_GO_FF)*M + params.δ.P*μᴾ*P + params.δ.D*μᴰ*D + params.λₚₒ*POC + (1 - params.γ.M)*Rᵤₚᴹ - Remin - Denit - ϕ₁ᴰᴼᶜ - ϕ₂ᴰᴼᶜ - ϕ₃ᴰᴼᶜ
end