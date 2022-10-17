#has to be discrete as needs access to field derivitatives
#I think the plankton sinking terms should be multiplied by the shear in order to conserve mass (not reflected in eqs 37 and 40) so going todo that
function POC_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #Microzooplankton feeding waste
    prey =  (; P, D, POC)

    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    Σgᶻ = gᶻₚ + gᶻ_D + gᶻₚₒ

    f_Z = params.b.Z^T

    z_waste = params.σ.Z*Σgᶻ*Z + params.r.Z*f_Z*Z^2/(params.Kₘ+Z) + params.m.Z*f_Z*Z^2

    #Phytoplankton
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    Kₙₒ₃, Kₙₕ₄ = K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄) 
    p_waste = (1 - 0.5*R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, Kₙₕ₄))*(params.m.P*P^2/(P+params.Kₘ)+sh*params.wᴾ*P^2)

    #aggregation
    ϕ₁ᴰᴼᶜ = sh*(params.a₁*DOC + params.a₂*POC)*DOC 
    ϕ₃ᴰᴼᶜ = (params.a₄*POC + params.a₅*DOC)*DOC

    ϕ = sh*params.a₆*POC^2 + sh*params.a₇*POC*GOC + params.a₈*POC*GOC + params.a₉*POC^2

    #mezozooplankton grazing
    prey = (; P, D, Z, POC)
    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ_FF = params.g_FF*params.b.M^T*params.wₚₒ*POC

    #microzooplankton grazing
    prey =  (; P, D, POC)

    L_day = params.L_day(t)
    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    return z_waste + 0.5*params.m.D*D^2/(D+params.Kₘ) + p_waste + params.λₚₒ*GOC + ϕ₁ᴰᴼᶜ + ϕ₃ᴰᴼᶜ - (gᴹₚₒ + gᴹₚₒ_FF)*M - gᶻₚₒ*Z - params.λₚₒ*POC -  ϕ# - params.wₚₒ*∂zᶜᶜᶜ(i, j, k, grid, model_fields.POC)
end

function GOC_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #Macrozooplankton feeding waste
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

    f_M = params.b.M^T

    w_GOC = (params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)
    gᴹₚₒ_FF = params.g_FF*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_FF*f_M*w_GOC*GOC

    m_waste = params.σ.M*(Σgᴹ + gᴹₚₒ_FF + gᴹ_GO_FF)*M +params.r.M*f_M*M^2/(M+params.Kₘ)

    #Upper trophic mesozooplankton excretion
    if params.upper_trophic_feeding
        Pᵤₚᴹ = params.σ.M*(1/(1-params.eₘₐₓ.M))*params.m.M*f_M*M^2
    else
        Pᵤₚᴹ = 0.0
    end

    #Phytoplankton
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    Kₙₒ₃, Kₙₕ₄ = K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄) 
    p_waste = (1 - 0.5*R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, Kₙₕ₄))*(params.m.P*P^2/(P+params.Kₘ)+sh*params.wᴾ*P^2)

    #Diatom waste
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Kₙₒ₃, Kₙₕ₄ = K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ)
    Lₙₒ₃ᴰ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴰ = Lₙₒ₃ᴰ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)  
    L_Feᴰ = L_Fe(D, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)
    Lₛᵢᴰ = L_mondo(Si, params.Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2/(params.Kₛᵢ^2+Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴰ = min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ)
    d_waste = 0.5*params.m.D*D^2/(D+params.Kₘ) + sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D^2 #I think there is a typo on line 4 of eq 40, looks correct on line 6/7 of eq 49

    #aggregation
    ϕ₂ᴰᴼᶜ = sh*params.a₃*GOC*DOC

    ϕ = sh*params.a₆*POC^2 + sh*params.a₇*POC*GOC + params.a₈*POC*GOC + params.a₉*POC^2

    return m_waste + Pᵤₚᴹ + p_waste + d_waste + ϕ + ϕ₂ᴰᴼᶜ  - gᴹ_GO_FF*M - params.λₚₒ*GOC# - w_GOC*∂zᶜᶜᶜ(i, j, k, grid, model_fields.GOC)
end

function Feᴾᴼ_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #Microzooplankton feeding waste
    prey =  (; P, D, POC)

    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚ = g(Fᵢ(params.p.Z.P, params.Jₜₕᵣ.Z, P), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 
    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    Σθ_Fegᶻ = gᶻₚ*Feᴾ/(P+eps(0.0))+ gᶻ_D*Feᴰ/(D+eps(0.0))+ gᶻₚₒ*Feᴾᴼ/(POC+eps(0.0)) 

    f_Z = params.b.Z^T

    z_waste = params.σ.Z*Σθ_Fegᶻ*Z + params.θᶠᵉ.Z*(params.r.Z*f_Z*Z^2/(params.Kₘ+Z) + params.m.Z*f_Z*Z^2)

    #Phytoplankton
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    Kₙₒ₃, Kₙₕ₄ = K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄) 
    p_waste = (1 - 0.5*R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, Kₙₕ₄))*(params.m.P*P^2/(P+params.Kₘ)+sh*params.wᴾ*P^2)

    #mezozooplankton grazing
    prey = (; P, D, Z, POC)
    gₘᴹ = params.gₘₐₓ⁰.M*params.b.M^T

    Fₗᵢₘᴹ = Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, prey)
    Fᴹ = F(params.p.M, params.Jₜₕᵣ.M, prey)
    food_totalᴹ = food_total(params.p.M, prey)

    gᴹₚₒ = g(params.p.M.POC*(POC-params.Jₜₕᵣ.M), gₘᴹ, Fₗᵢₘᴹ, Fᴹ, params.K_G.M, food_totalᴹ) 
    gᴹₚₒ_FF = params.g_FF*params.b.M^T*params.wₚₒ*POC

    #microzooplankton grazing
    prey =  (; P, D, POC)

    L_day = params.L_day(t)
    gₘᶻ = params.gₘₐₓ⁰.Z*params.b.Z^T

    Fₗᵢₘᶻ = Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, prey)
    Fᶻ = F(params.p.Z, params.Jₜₕᵣ.Z, prey)
    food_totalᶻ = food_total(params.p.Z, prey)

    gᶻₚₒ = g(params.p.Z.POC*(POC-params.Jₜₕᵣ.Z), gₘᶻ, Fₗᵢₘᶻ, Fᶻ, params.K_G.Z, food_totalᶻ) 

    #aggregation
    ϕ = sh*params.a₆*POC^2 + sh*params.a₇*POC*GOC + params.a₈*POC*GOC + params.a₉*POC^2

    #Bacteria
    bFe = Fe

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)

    Bactfe = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))*BactfeR(params.μₘₐₓ⁰*params.bₚ^T, Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.Kₙₒ₃.Bact, params.Kₙₕ₄.Bact), bFe, params.K_Fe.Bact, params.θₘₐₓᶠᵉᴮᵃᶜᵗ)
    return (z_waste + 
        (Feᴰ/(D+eps(0.0)))*0.5*params.m.D*D^2/(D+params.Kₘ) + 
        p_waste + 
        params.λ_Fe*POC*Feᶠ(Fe, DOC, T) + 
        Cgfe1(DOC, POC, Fe, T, sh, params.a₁, params.a₂, params.a₄, params.a₅) - 
        params.λₚₒ*Feᴾᴼ + 
        params.κᶠᵉᴾᴼ_Bact*Bactfe - 
        ((gᴹₚₒ + gᴹₚₒ_FF)*M + gᶻₚₒ*Z)*Feᴾᴼ/(POC+eps(0.0)) - 
        (Feᴾᴼ/(POC+eps(0.0)) )*ϕ# - 
        #params.wₚₒ*∂zᶜᶜᶜ(i, j, k, grid, model_fields.Feᴾᴼ))
    )
end


function Feᴳᴼ_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)

    #Macrozooplankton feeding waste
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
    Σθ_Fegᴹ = gᴹₚ*Feᴾ/(P+eps(0.0))+ gᴹ_D*Feᴰ/(D+eps(0.0))+ gᴹₚₒ*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_Z*params.θᶠᵉ.Z
 
    f_M = params.b.M^T
 
    gᴹₚₒ_FF = params.g_FF*f_M*params.wₚₒ*POC
    gᴹ_GO_FF = params.g_FF*f_M*(params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)*GOC

    #Upper trophic mesozooplankton excretion
    if params.upper_trophic_feeding
        Pᵤₚᴹ = params.σ.M*(1/(1-params.eₘₐₓ.M))*params.m.M*f_M*M^2
    else
        Pᵤₚᴹ = 0.0
    end
 
    m_waste = params.σ.M*(Σθ_Fegᴹ + gᴹₚₒ_FF*Feᴾᴼ/(POC+eps(0.0)) + gᴹ_GO_FF*Feᴳᴼ/(GOC+eps(0.0)) )*M + params.θᶠᵉ.M*(params.r.M*f_M*M^2/(M+params.Kₘ) + Pᵤₚᴹ)

    #Phytoplankton
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    Kₙₒ₃, Kₙₕ₄ = K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ)
    L_NO₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴾ = L_NO₃ᴾ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄) 
    p_waste = (1 - 0.5*R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, params.r_CaCO₃, Lₙᴾ, Kₙₕ₄))*(params.m.P*P^2/(P+params.Kₘ)+sh*params.wᴾ*P^2)*Feᴾ/(P+eps(0.0))

    #Diatom waste
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D

    Lₚₒ₄ᴰ = L_mondo(PO₄, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ))
    Kₙₒ₃, Kₙₕ₄ = K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ)
    Lₙₒ₃ᴰ = L_NO₃(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)
    Lₙᴰ = Lₙₒ₃ᴰ + L_NH₄(NO₃, NH₄, Kₙₒ₃, Kₙₕ₄)  
    L_Feᴰ = L_Fe(D, Chlᴰ, Feᴰ, params.θᶠᵉₒₚₜ.D, Lₙᴰ, Lₙₒ₃ᴰ)
    Lₛᵢᴰ = L_mondo(Si, params.Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2/(params.Kₛᵢ^2+Si̅^2)) # eq 11b, 12. Need a callback that records Sī

    Lₗᵢₘᴰ = min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ)
    d_waste = (0.5*params.m.D*D^2/(D+params.Kₘ) + sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D^2)*Feᴰ/(D+eps(0.0))

    #sinking
    w_GOC = (params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)

    #grazing by zooplankton
    gᴹ_GO_FF = params.g_FF*f_M*w_GOC*GOC*Feᴳᴼ/(GOC+eps(0.0)) 

    #aggregation
    ϕ = sh*params.a₆*POC^2 + sh*params.a₇*POC*GOC + params.a₈*POC*GOC + params.a₉*POC^2

    #Bacteria
    bFe = Fe

    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.P, params.Kₙₒ₃ᵐⁱⁿ.P, params.Kₙₕ₄ᵐⁱⁿ.P, params.Sᵣₐₜ.P, params.Pₘₐₓ.P
    μᴾ, Lₗᵢₘᴾ = μᵢ(P, Chlᴾ, Feᴾ, Siᴰ, NO₃, NH₄, PO₄, 1.0, PARᴾ, T, zₘₓₗ, zₑᵤ, K(P, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), K(P, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Pₘₐₓ), params.θᶠᵉₒₚₜ.P, 0.0, 0.0, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.P, L_day, params.α.P)
    
    Bactfe = Bact(x, y, z, model_fields.Z, model_fields.M, max(zₘₓₗ, zₑᵤ))*BactfeR(params.μₘₐₓ⁰*params.bₚ^T, Lₗᵢₘᴮᵃᶜᵗ(bFe, PO₄, NO₃, NH₄, params.K_Fe.Bact, params.Kₚₒ₄.Bact, params.Kₙₒ₃.Bact, params.Kₙₕ₄.Bact), bFe, params.K_Fe.Bact, params.θₘₐₓᶠᵉᴮᵃᶜᵗ)

    return (m_waste + p_waste + d_waste + gᴹ_GO_FF*M +
        ϕ*Feᴾᴼ/(POC+eps(0.0)) + 
        params.λ_Fe*GOC*Feᶠ(Fe, DOC, T) +
        Cgfe2(DOC, GOC, Fe, T, sh, params.a₃) +
        params.κᶠᵉᴳᴼ_Bact*Bactfe -
        params.λₚₒ*Feᴳᴼ# - 
        #w_GOC*∂zᶜᶜᶜ(i, j, k, grid, model_fields.Feᴳᴼ))
    )
end

function Siᴾ_forcing(i, j, k, grid, clock, model_fields, params)
    P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₘₓₗ, zₑᵤ, Si̅, D_dust = get_local_value.(i, j, k, values(model_fields[(:P, :D, :Chlᴾ, :Chlᴰ, :Feᴾ, :Feᴰ, :Siᴰ, :Z, :M, :DOC, :POC, :GOC, :Feᴾᴼ, :Feᴳᴼ, :Siᴾ, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :PAR¹, :PAR², :PAR³, :T, :S, :zₘₓₗ, :zₑᵤ, :Si̅, :D_dust)]))
    x, y, z, t = grid.xᶜᵃᵃ[i], grid.yᵃᶜᵃ[j], grid.zᵃᵃᶜ[k], clock.time
    PARᴾ, PARᴰ, PAR = PAR_components(PAR¹, PAR², PAR³, params.β₁, params.β₂, params.β₃)
    sh = get_sh(z, zₘₓₗ, params.shₘₓₗ, params.shₛᵤ)
    L_day = params.L_day(t)

    #diatom grazing
    z_food_total = food_total(params.p.Z, (; P, D, POC)) 
    m_food_total = food_total(params.p.M, (; P, D, POC, Z)) 
    
    gᶻ_D = g(Fᵢ(params.p.Z.D, params.Jₜₕᵣ.Z, D), params.gₘₐₓ⁰.Z*params.b.Z^T, Fₗᵢₘ(params.p.Z, params.Jₜₕᵣ.Z, params.Fₜₕᵣ.Z, (; P, D, POC)), F(params.p.Z, params.Jₜₕᵣ.Z, (; P, D, POC)), params.K_G.Z, z_food_total)
    gᴹ_D = g(Fᵢ(params.p.M.D, params.Jₜₕᵣ.M, D), params.gₘₐₓ⁰.M*params.b.M^T, Fₗᵢₘ(params.p.M, params.Jₜₕᵣ.M, params.Fₜₕᵣ.M, (; P, D, POC, Z)), F(params.p.M, params.Jₜₕᵣ.M, (; P, D, POC, Z)), params.K_G.M, m_food_total)
    
    #diatom mortality
    Kₚₒ₄ᵐⁱⁿ, Kₙₒ₃ᵐⁱⁿ, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ = params.Kₚₒ₄ᵐⁱⁿ.D, params.Kₙₒ₃ᵐⁱⁿ.D, params.Kₙₕ₄ᵐⁱⁿ.D, params.Sᵣₐₜ.D, params.Pₘₐₓ.D
    μᴰ, Lₗᵢₘᴰ  = μᵢ(D, Chlᴰ, Feᴰ, Siᴰ, NO₃, NH₄, PO₄, Si, PARᴰ, T, zₘₓₗ, zₑᵤ, K(D, Kₚₒ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₒ₃ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), K(D, Kₙₕ₄ᵐⁱⁿ, Sᵣₐₜ, Dₘₐₓ), params.θᶠᵉₒₚₜ.D, params.Kₛᵢᴰᵐⁱⁿ, Si̅, params.Kₛᵢ, params.μₘₐₓ⁰, params.bₚ, params.t_dark.D, L_day, params.α.D)

    d_waste = (params.m.D*L_mondo(D, params.Kₘ) + sh*(params.wᴾ+ params.wₘₐₓᴰ*(1-Lₗᵢₘᴰ))*D)*Siᴰ #I think there is a mistake in the second term of line 2 eq 51 and the \theta  shouldn't be there, sometimes they write X^2*θʸˣ and sometimes X*Yˣ even though they're the same thing

    #sinking 
    w_GOC = (params.w_GOᵐⁱⁿ + (200-params.w_GOᵐⁱⁿ)*max(0, -z - max(zₑᵤ, zₘₓₗ))/5000)

    #dissolusion
    zₘₐₓ = max(zₑᵤ, zₘₓₗ)
    χₗₐᵦ = params.χₗₐᵦ⁰*ifelse(z<=zₘₐₓ, 1.0, exp(-(params.λₚₛᵢˡᵃᵇ - params.λₚₛᵢʳᵉᶠ)*((z-zₘₐₓ)/w_GOC)))
    λₚₛᵢ = χₗₐᵦ*params.λₚₛᵢˡᵃᵇ+(1-χₗₐᵦ)*params.λₚₛᵢʳᵉᶠ

    Siₑ = 10^(6.44-968/(T+273.15))
    Siₛₐₜ = (Siₑ - Si)/Siₑ
    λₚₛᵢ *= *(0.225*(1+T/15)*Siₛₐₜ + 0.775*((1+T/400)^4*Siₛₐₜ)^9)

    #Dissₛₛ is not defined anywhere so I think it must be a typo and not be anythng?
    #removing it is more concistent with the formulation of the other particle properties

    return (Siᴰ/(D+eps(0.0)))*(gᶻ_D*Z + gᴹ_D*M) + d_waste - λₚₛᵢ*Siᴾ# - w_GOC*∂zᶜᶜᶜ(i, j, k, grid, model_fields.CaCO₃) 
end