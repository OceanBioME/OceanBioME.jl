@inline function Bact(x, y, z, Z, M, zₘₐₓ)
    #Do macroalgaes not contribute to the bacteria pool? So need to add some proxy from their biomass
    Zᵢⱼₖ = interpolate(Z, x, y, z)
    Mᵢⱼₖ = interpolate(M, x, y, z)
    return ifelse(z<=zₘₐₓ, min(0.7*(Zᵢⱼₖ + 2*Mᵢⱼₖ), 4), Bact(x, y, zₘₐₓ, k, Z, M, zₘₐₓ)*(zₘₐₓ/z)^0.683)
end

#since Bact needs access to the whole column of Z and M need to use discrete forcings anything that depends on it
function DOC_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾᴼ, Siᴳᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, PARᴾ, PARᴰ, T, zₘₓₗ, zₑᵤ, ϕ = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾᴼ, :Siᴳᴼ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :O₂, :PARᴾ, :PARᴰ, :T, :zₘₓₗ, :zₑᵤ, :ϕ)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock

    L_day = params.L_day(t)
    #Plankton quality
    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Inf, PARᴾ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P, params.K_NO₃.P, params.K_NH₄.P, params.θᶠᵉₒₚₜ.P, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D, params.K_NO₃.D, params.K_NH₄.D, params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #Microzooplankton efficiency
    prey =  (; P, D, POC)

    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(params.p.Z.P*(P-params.Jₜₕᵣ.Z.P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(params.p.Z.D*(D-params.Jₜₕᵣ.Z.D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z.POC), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    Σgᶻ = gᶻₚ + gᶻ_D + gᶻₚₒ

    gᶻ_Fe = gᶻₚ*Feᴾ/P + gᶻ_D*Feᴰ/D + gᶻₚₒ*Feᴾᴼ/POC
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*(gᶻₚ + gᶻ_D + gᶻₚₒ)), gᶻ_Fe/(params.θᶠᵉ.Z*(gᶻₚ + gᶻ_D + gᶻₚₒ)))
    eᶻ = eₙᶻ*min(params.eₘₐₓᶻ, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*(gᶻₚ + gᶻ_D + gᶻₚₒ)))

    #Mezozooplankton efficiency
    prey =  (; P, D, POC, Z)

    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚ = g(params.p.M.P*(P-params.Jₜₕᵣ.M.P), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_D = g(params.p.M.D*(D-params.Jₜₕᵣ.M.D), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M.POC), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_Z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.M.Z), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ)

    f_M = params.b.M^T
    
    Σgᴹ = gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹ_Z

    gᴹ_Fe = gᴹₚ*Feᴾ/P + gᴹ_D*Feᴰ/D + gᴹₚₒ*Feᴾᴼ/POC + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z)*params.θᴺᶜ

    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ, gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ)))
    eᴹ = eₙᴹ*min(params.eₘₐₓᴹ, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ))

    gᴹ_GO_FF = params.g_ff*f_M*(params._GOᵐⁱⁿ + (200-params._GOᵐⁱⁿ)*max(0, z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    #plankton growth
    μᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Inf, PARᴾ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P, params.K_NO₃.P, params.K_NH₄.P, params.θᶠᵉₒₚₜ.P, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    μᴰ = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D, params.K_NO₃.D, params.K_NH₄.D, params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    #Upper trophic mesozooplankton feeding
    if params.upper_trophic_feeding
        Rᵤₚᴹ = (1 - params.eₘₐₓᴹ- params.σᴹ)*(1/(1-params.eₘₐₓᴹ))*params.m.M*f_M*M^2
    else
        Rᵤₚᴹ = 0.0
    end

    #Bacterial action
    Bactᵢⱼₖ = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))

    bFe = Fe
    L_Feᴮᵃᶜᵗ = L_mondo(Fe, params.K_Fe.Bact)
    Lₚₒ₄ᴮᵃᶜᵗ = L_mondo(PO₄, params.Kₚₒ₄ᴮᵃᶜᵗ)
    L_NO₃ᴮᵃᶜᵗ = L_NO₃(NO₃, NH₄, params.K_NO₃.Bact, params.K_NH₄.Bact)
    Lₙᴮᵃᶜᵗ = L_NO₃ᴮᵃᶜᵗ + L_NH₄(NO₃, NH₄, params.K_NO₃.Bact, params.K_NH₄.Bact) 
    Lₗᵢₘᴮᵃᶜᵗ = min(Lₙᴮᵃᶜᵗ, Lₚₒ₄ᴮᵃᶜᵗ, L_Feᴮᵃᶜᵗ)#eq 34c, assuming typo of Lₙₕ₄ vs Lₙ because why else is it defined?
    
    Lᴮᵃᶜᵗ = Lₗᵢₘᴮᵃᶜᵗ*L_mondo(DOC, params.K_DOC)
    μ_Bact = params.λ_DOC*params.bₚ*Bact*DOC*Lᴮᵃᶜᵗ/params.Bactᵣₑ #eq 33a, b say Lₗᵢₘᴮᵃᶜᵗ but think this must be Lᴮᵃᶜᵗ

    Remin = min(O₂/params.O₂ᵘᵗ, μ_Bact*(1-ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)))
    Denit = min(NO₃/params.rₙₒ₃⋆, μ_Bact*ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²))

    #aggregation
    sh = if(zₘₓₗ<z, params.shₘₓₗ, params.shₛᵤ)
    
    ϕ₁ᴰᴼᶜ = sh*(params.a₁*DOC + params.a₂*POC)*DOC 
    ϕ₂ᴰᴼᶜ = sh*params.a₃*GOC*DOC
    ϕ₃ᴰᴼᶜ = (params.a₄*POC + params.a₅*DOC)*DOC
    return (1-params.γ.Z)*(1 - eᶻ - params.σᶻ)*Σgᶻ*Z + (1-γ.M)*(1 - eᴹ - params.σᴹ)*(Σgᴹ + gᴹ_GO_FF)*M + params.δ.P*μᴾ*P + params.δ.D*μᴰ*D + λₚₒ⋆*POC + (1 - params.γ.M)*Rᵤₚᴹ - Remin - Denit - ϕ₁ᴰᴼᶜ - ϕ₂ᴰᴼᶜ - ϕ₃ᴰᴼᶜ
end