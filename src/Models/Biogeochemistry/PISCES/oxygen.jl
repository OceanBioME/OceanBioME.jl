
function O₂_forcing(i, j, k, grid, clock, model_fields, params) #needs to be discrete because of \bar{PAR}
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾᴼ, Siᴳᴼ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, ϕ = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾᴼ, :Siᴳᴼ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :ϕ)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)
    L_day = params.L_day(t)

    #phytoplankton growth
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Inf, PARᴾ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P, params.K_NO₃.P, params.K_NH₄.P, params.θᶠᵉₒₚₜ.P, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D, params.K_NO₃.D, params.K_NH₄.D, params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)
    
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P)
    L_NH₄ᴾ =  + L_NH₄(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P)
    μₙₕ₄ᴾ = μᴾ * L_NH₄ᴾ / (L_NO₃ᴾ + L_NH₄ᴾ)
    μₙₒ₃ᴾ = μᴾ * L_NO₃ᴾ / (L_NO₃ᴾ + L_NH₄ᴾ)

    L_NO₃ᴰ = L_NO₃(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D)
    L_NH₄ᴰ =  + L_NH₄(NO₃, NH₄, params.K_NO₃.D, params.K_NH₄.D)
    μₙₕ₄ᴰ = μᴰ * L_NH₄ᴰ / (L_NO₃ᴰ + L_NH₄ᴰ)
    μₙₒ₃ᴰ = μᴰ * L_NO₃ᴰ / (L_NO₃ᴰ + L_NH₄ᴰ)

    #microzooplankton grazing
    prey =  (; P, D, POC)
    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(params.p.Z.P*(P-params.Jₜₕᵣ.Z.P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(params.p.Z.D*(D-params.Jₜₕᵣ.Z.D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z.POC), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    Σgᶻ = gᶻₚ + gᶻ_D + gᶻₚₒ

    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Inf, PARᴾ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P, params.K_NO₃.P, params.K_NH₄.P, params.θᶠᵉₒₚₜ.P, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D, params.K_NO₃.D, params.K_NH₄.D, params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    gᶻ_Fe = gᶻₚ*Feᴾ/P + gᶻ_D*Feᴰ/D + gᶻₚₒ*Feᴾᴼ/POC
    gᶻₙ = gᶻₚ*θᴺᴾ + gᶻ_D*θᴺᴰ + gᶻₚₒ*params.θᴺᶜ

    eₙᶻ = min(1, gᶻₙ/(params.θᴺᶜ*Σgᶻ), gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ))
    eᶻ = eₙᶻ*min(params.eₘₐₓᶻ, (1 - params.σ.Z)*gᶻ_Fe/(params.θᶠᵉ.Z*Σgᶻ))

    #Macrozooplankton grazing
    prey =  (; P, D, POC, Z)

    L_day = params.L_day(t)
    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚ = g(params.p.M.P*(P-params.Jₜₕᵣ.M.P), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_D = g(params.p.M.D*(D-params.Jₜₕᵣ.M.D), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M.POC), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹ_Z = g(params.p.M.Z*(Z - params.Jₜₕᵣ.M.Z), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ)
    Σgᴹ = gᴹₚ + gᴹ_D + gᴹₚₒ + gᴹ_Z

    θᴺᴾ = θᴺᴵ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, Inf, PARᴾ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P, params.K_NO₃.P, params.K_NH₄.P, params.θᶠᵉₒₚₜ.P, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    θᴺᴰ = θᴺᴵ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, params.Kᵐⁱⁿₚₒ₄.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D, params.K_NO₃.D, params.K_NH₄.D, params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, params.Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    f_M = params.b.M^T

    gᴹₚₒ_FF = params.g_ff*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_ff*f_M*(params.w_GOᵐⁱⁿ + (200-params._GOᵐⁱⁿ)*max(0, z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    gᴹ_Fe = gᴹₚ*Feᴾ/P + gᴹ_D*Feᴰ/D + gᴹₚₒ*Feᴾᴼ/POC + gᴹ_Z*params.θᶠᵉ.Z
    gᴹₙ = gᴹₚ*θᴺᴾ + gᴹ_D*θᴺᴰ + (gᴹₚₒ + gᴹ_Z)*params.θᴺᶜ

    eₙᴹ = min(1, gᴹₙ/(params.θᴺᶜ*Σgᴹ), gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ))
    eᴹ = eₙᴹ*min(params.eₘₐₓᴹ, (1 - params.σ.M)*gᴹ_Fe/(params.θᶠᵉ.M*Σgᴹ))

    #Upper trophic feeding waste
    if params.upper_trophic_feeding
        Rᵤₚᴹ = (1 - params.σ.M - params.eₘₐₓᴹ)*(1/(1-params.eₘₐₓᴹ))*params.m.M*f_M*M^2
    else
        Rᵤₚᴹ = 0.0
    end

    #Bacteria
    ΔO₂ᵢⱼₖ = ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)

    Lᴮᵃᶜᵗ = Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.K_NO₃.Bact, params.K_NH₄.bact)*L_mondo(DOC, params.K_DOC)
    μ_Bact = params.λ_DOC*params.bₚ*Bact*DOC*Lᴮᵃᶜᵗ/params.Bactᵣₑ #eq 33a, b say Lₗᵢₘᴮᵃᶜᵗ but think this must be Lᴮᵃᶜᵗ

    Remin = min(O₂/params.O₂ᵘᵗ, μ_Bact*(1-ΔO₂(O₂, params.O₂ᵐⁱⁿ¹, params.O₂ᵐⁱⁿ²)))
    Denit = min(NO₃/params.rₙₒ₃, μ_Bact*ΔO₂ᵢⱼₖ)

    kₘₓₗ = round(Int, fractional_z_index(zₘₓₗ, Center(), grid))+1
    PAR̄ = sum(PAR[i, j, kₘₓₗ:grid.Nz])/zₘₓₗ
    Nitrif =  params.λ.NH₄*NH₄*(1-ΔO₂ᵢⱼₖ)/(1+PAR̄)

    μₚ = params.μₘₐₓ⁰*params.bₚ^T
    bFe = Fe

    Lₙᴾ = L_NO₃(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P) + L_NH₄(NO₃, NH₄, params.K_NO₃.P, params.K_NH₄.P)
    Lₙᴰᶻ = ifelse(Lₙᴾ >= 0.8, 0.01, 1 - Lₙᴾ)
    N_fix = params.N_fixᵐ*max(0, μₚ - 2.15)*Lₙᴰᶻ*min(L_mondo(bFe, params.K_Feᴰᶻ), L_mondo(PO₄, params.Kᵐⁱⁿₚₒ₄.P)*(1 - exp(-PAR/params.E_fix))

    return (
        params.O₂ᵘᵗ*(μₙₕ₄ᴾ*P + μₙₕ₄ᴰ*D)
        + (params.O₂ᵘᵗ + params.O₂ⁿⁱᵗ)*(μₙₒ₃ᴾ*P + μₙₒ₃ᴰ*D)
        - params.O₂ᵘᵗ*params.γ.Z*(1 - eᶻ - params.σᶻ)*Σgᶻ*Z 
        - params.O₂ᵘᵗ*params.γ.M*(1 - eᴹ - params.σ.M)*(Σgᴹ + gᴹₚₒ_FF + gᴹ_GO_FF)*M #must be a typo on line 4 of eq 81 since Rᵤₚᴹ isn't a per M value and shouldn't get multipled by θ and γ twice? 
        - params.O₂ᵘᵗ*params.γ.M*Rᵤₚᴹ 
        - params.O₂ᵘᵗ*(Remin + Denit)
        + params.O₂ᵘᵗ*N_fix
    ) 
end