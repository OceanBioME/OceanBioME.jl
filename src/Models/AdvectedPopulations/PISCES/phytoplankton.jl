
#PISCES is an implementation of PISCES-v2 described in the 2015 Aumont paper.
#Where possible conventions have been kept in line with the original paper, including notation and equation groupings. 
#Some changes have been made due to queries with the original paper, but changes and justification for these have been noted.

#This document contains equations for:
    #nutrient_quotas
    #concentration_limitation for determining Monod nutrient limitation terms

    #shear_rate
    #day_length
    #day_dependent_growth_rate (eq3a), depth_dependent_growth_rate (eq3d), t_dark (eq3b)
    #Half saturation constants
    #P_nutrient_limitation (eq6a), D_nutrient_limitation (eq11a)

    #minimum_iron_quota

    #P_PAR
    #phytoplankton_iron_biomass_growth_rate
    #phytoplankton_growth_rate
    #nutrient limitation
    #variation_in_SiC_ratio
    #D_quadratic_mortality
    #Forcing equations

@inline nutrient_quota(I, J) = ifelse(J == 0, 0, I/(J + eps(0.0))) # this shouldn't be needed as I should be zero if J is zero

@inline concentration_limitation(I, J) = I/(I + J + eps(0.0))

#Expresses growth rate with dependency on day length
@inline day_dependent_growth_rate(L_day) = 1.5*concentration_limitation(L_day, 0.5)  #eq 3a

#Mean time phytoplankton can spend in unlit part of mixed layer.
@inline function t_dark(zₘₓₗ, zₑᵤ, bgc)
    κᵥₑᵣₜ = bgc.vertical_diffusivity  
    return (max(0, abs(zₘₓₗ)-abs(zₑᵤ))^2)/(κᵥₑᵣₜ[0,0,0] + eps(0.0)) #eq 3b,c, do not divide by 86400. Vertical diffusivity 
end
@inline depth_dependent_growth_rate(zₘₓₗ, zₑᵤ, t_darkᴵ, bgc) = 1 - concentration_limitation(t_dark(zₘₓₗ, zₑᵤ, bgc), t_darkᴵ) #eq 3d

#The minimum iron quota is the sum of the three demands for iron in phytoplankton (photosynthesis, respiration, nitrate reduction)
@inline minimum_iron_quota(I, Iᶜʰˡ, Lₙᴵ, Lₙₒ₃ᴵ) = 0.0016/(55.85) * nutrient_quota(Iᶜʰˡ, I) + 1.21e-5*14*Lₙᴵ/(55.85*7.625)*1.5+1.15e-4*14*Lₙₒ₃ᴵ/(55.85*7.625) #eq 20 -> Lₙ could be meant to be L_NH₄?

#Different size classes of phytoplankton, have different half-saturation constants due to varying surface area to volume ratios. 
#Generally increased biomass corresponds to larger size classes. 
#Half saturation constants have biomass dependency.
@inline I₁(I, Iₘₐₓ) = min(I, Iₘₐₓ) #eq 7a
@inline I₂(I, Iₘₐₓ) = max(0, I - Iₘₐₓ) #eq 7b
@inline nutrient_half_saturation_const(Kᵢᴶᵐⁱⁿ, J₁, J₂, Sᵣₐₜᴶ) = Kᵢᴶᵐⁱⁿ* (J₁ + Sᵣₐₜᴶ* J₂)/(J₁ + J₂ + eps(0.0)) #eq 7c

#Light absorption by phytoplankton. Visible light split into 3 wavebands, where light absorption of each waveband controlled by coefficient.
@inline function P_PAR(PAR₁, PAR₂, PAR₃, bgc)
    β₁ᴾ = bgc.absorption_in_the_blue_part_of_light.P
    β₂ᴾ = bgc.absorption_in_the_green_part_of_light.P
    β₃ᴾ = bgc.absorption_in_the_red_part_of_light.P

    return β₁ᴾ*PAR₁ + β₂ᴾ*PAR₂ + β₃ᴾ*PAR₃
end

@inline function D_PAR(PAR₁, PAR₂, PAR₃, bgc)
    β₁ᴰ = bgc.absorption_in_the_blue_part_of_light.D
    β₂ᴰ = bgc.absorption_in_the_green_part_of_light.D
    β₃ᴰ = bgc.absorption_in_the_red_part_of_light.D

    return β₁ᴰ*PAR₁ + β₂ᴰ*PAR₂ + β₃ᴰ*PAR₃
end

#The growth rate of the iron biomass of phytoplankton.
@inline function phytoplankton_iron_biomass_growth_rate(I, Iᶠᵉ, θₘₐₓᶠᵉᴵ, Sᵣₐₜᴵ, K_Feᴵᶠᵉᵐⁱⁿ, Iₘₐₓ, L_Feᴵ, bFe, T, bgc) 
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    bₚ = bgc.temperature_sensitivity_of_growth

    μₚ = μ⁰ₘₐₓ*(bₚ^T) #4b

    I₂ = max(0, I - Iₘₐₓ) #18c
    I₁ = I - I₂ #18c

    K_Feᴵᶠᵉ = K_Feᴵᶠᵉᵐⁱⁿ*(I₁ + Sᵣₐₜᴵ*I₂)/(I₁+I₂+eps(0.0)) #18b

    Lₗᵢₘ₁ᴵᶠᵉ = concentration_limitation(bFe, K_Feᴵᶠᵉ) #18a
    #Lₗᵢₘ₂ᴵᶠᵉ = (4 - 4.5*L_Feᴵ)/(L_Feᴵ + 0.5) #Formulation given in paper does not vary between 1 and 4 as claimed
    Lₗᵢₘ₂ᴵᶠᵉ = (4 - 2*L_Feᴵ)/(L_Feᴵ + 1) #19, reformulation is between given bounds

    return θₘₐₓᶠᵉᴵ*Lₗᵢₘ₁ᴵᶠᵉ*Lₗᵢₘ₂ᴵᶠᵉ*(1 - (nutrient_quota(Iᶠᵉ, I))/(θₘₐₓᶠᵉᴵ + eps(0.0)))/(1.05 - (nutrient_quota(Iᶠᵉ, I))/(θₘₐₓᶠᵉᴵ + eps(0.0)))*μₚ #eq17
end

#Growth rates of phytoplankton depend on limiting nutrients, day length, and light absorbtion of phytoplankton.
@inline function phytoplankton_growth_rate(I, Iᶜʰˡ, PARᴵ, L_day, T, αᴵ, Lₗᵢₘᴵ, zₘₓₗ, zₑᵤ, t_darkᴵ, bgc)
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    bₚ = bgc.temperature_sensitivity_of_growth

    μₚ = μ⁰ₘₐₓ*(bₚ^T)    #eq 4b      

    return μₚ * day_dependent_growth_rate(L_day) * depth_dependent_growth_rate(zₘₓₗ, zₑᵤ, t_darkᴵ, bgc) * (1-exp(-αᴵ*(nutrient_quota(Iᶜʰˡ,I))*PARᴵ/(L_day*μₚ*Lₗᵢₘᴵ + eps(0.0)))) * Lₗᵢₘᴵ #eq2b 
end


#Nutrient limitation terms.
#Nutrient and phosphate limitations are based on Monod parametrisations, iron on quota parametrisations.
@inline ammonium_limitation(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₒ₃ᴵ*NH₄/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄ + eps(0.0)) #eq 6d
@inline nitrate_limitation(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₕ₄ᴵ*NO₃/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄ + eps(0.0)) #eq 6e
@inline iron_limitation(I, Iᶠᵉ, θₒₚₜᶠᵉᴵ, θₘᵢₙᶠᵉᴵ) = min(1, max(0, (nutrient_quota(Iᶠᵉ, I) - θₘᵢₙᶠᵉᴵ)/(θₒₚₜᶠᵉᴵ + eps(0.0)))) #eq 6f

#Determines individual nutrient limitation terms, and overall limiting nutrients.
@inline function P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)
    #Parameters
    θₒₚₜᶠᵉᵖ = bgc.optimal_iron_quota.P
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton.P
    Kₙₒ₃ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_nitrate.P
    Kₙₕ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium.P
    Kₚₒ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate.P 
    Pₘₐₓ = bgc.threshold_concentration_for_size_dependency.P

    #Half saturation constants
    P₁ = I₁(P, Pₘₐₓ)
    P₂ = I₂(P, Pₘₐₓ)
    Kₙₒ₃ᴾ = nutrient_half_saturation_const(Kₙₒ₃ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)
    Kₙₕ₄ᴾ = nutrient_half_saturation_const(Kₙₕ₄ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)
    Kₚₒ₄ᴾ = nutrient_half_saturation_const(Kₚₒ₄ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)

    #Nutrient limitation terms (for phosphate, ammonium, nitrate, iron)
    Lₚₒ₄ᴾ = concentration_limitation(PO₄, Kₚₒ₄ᴾ) #6b
    Lₙₕ₄ᴾ = ammonium_limitation(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙₒ₃ᴾ = nitrate_limitation(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙᴾ = Lₙₒ₃ᴾ + Lₙₕ₄ᴾ #6c
    θₘᵢₙᶠᵉᵖ = minimum_iron_quota(P, Pᶜʰˡ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ) #changed from Lₙᴾ to Lₙₕ₄ᴾ
    L_Feᴾ = iron_limitation(P, Pᶠᵉ, θₒₚₜᶠᵉᵖ, θₘᵢₙᶠᵉᵖ)

    return min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ), Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ #6a
end

#Same for Lₗᵢₘᴰ
@inline function D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    #Parameters
    θₒₚₜᶠᵉᴰ = bgc.optimal_iron_quota.D
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton.D
    Kₙₒ₃ᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_nitrate.D
    Kₙₕ₄ᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium.D
    Kₚₒ₄ᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate.D
    Kₛᵢᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_silicate
    Kₛᵢ = bgc.parameter_for_half_saturation_const
    Dₘₐₓ = bgc.threshold_concentration_for_size_dependency.D

    #Half saturation constants
    D₁ = I₁(D, Dₘₐₓ)
    D₂ = I₂(D, Dₘₐₓ)
    Kₙₒ₃ᴰ = nutrient_half_saturation_const(Kₙₒ₃ᴰᵐⁱⁿ, D₁, D₂, Sᵣₐₜᴰ)
    Kₙₕ₄ᴰ = nutrient_half_saturation_const(Kₙₕ₄ᴰᵐⁱⁿ, D₁, D₂, Sᵣₐₜᴰ)
    Kₚₒ₄ᴰ = nutrient_half_saturation_const(Kₚₒ₄ᴰᵐⁱⁿ, D₁, D₂, Sᵣₐₜᴰ)

    #Nutrient limitation terms (for phosphate, ammonium, nitrate, iron, silicate)
    Lₚₒ₄ᴰ = concentration_limitation(PO₄, Kₚₒ₄ᴰ) #6b
    Lₙₕ₄ᴰ = ammonium_limitation(NO₃, NH₄, Kₙₒ₃ᴰ, Kₙₕ₄ᴰ)
    Lₙₒ₃ᴰ = nitrate_limitation(NO₃, NH₄, Kₙₒ₃ᴰ, Kₙₕ₄ᴰ)
    Lₙᴰ = Lₙₒ₃ᴰ + Lₙₕ₄ᴰ         #6c
    θₘᵢₙᶠᵉᴰ = minimum_iron_quota(D, Dᶜʰˡ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ) #changed from n to NH₄
    L_Feᴰ = iron_limitation(D, Dᶠᵉ ,θₒₚₜᶠᵉᴰ, θₘᵢₙᶠᵉᴰ)
    Kₛᵢᴰ = Kₛᵢᴰᵐⁱⁿ + 7*Si̅^2 / (Kₛᵢ^2 + Si̅^2 + eps(0.0)) #12
    Lₛᵢᴰ = concentration_limitation(Si, Kₛᵢᴰ)    #11b

    return min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ), Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ,  L_Feᴰ, Lₛᵢᴰ   #11a
end


@inline function variation_in_SiC_ratio(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ, Si̅, bgc)
    θₘˢⁱᴰ = bgc.optimal_SiC_uptake_ratio_of_diatoms
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    Kₛᵢ¹ = bgc.parameter_for_SiC.one
    Kₛᵢ² = bgc.parameter_for_SiC.two
    bₚ = bgc.temperature_sensitivity_of_growth

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    
    μₚ = μ⁰ₘₐₓ*(bₚ^T)
    
    Lₗᵢₘ₁ᴰˢⁱ = concentration_limitation(Si, Kₛᵢ¹)    #23c
    Lₗᵢₘ₂ᴰˢⁱ = ifelse(ϕ < 0, (concentration_limitation((Si)^3, (Kₛᵢ²)^3)), 0)   #23d
    
    Fₗᵢₘ₁ᴰˢⁱ = min((μᴰ)/(μₚ*Lₗᵢₘᴰ + eps(0.0)), Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ)  #23a
    Fₗᵢₘ₂ᴰˢⁱ = min(1, 2.2*max(0, Lₗᵢₘ₁ᴰˢⁱ - 0.5)) #23b

    return θₘˢⁱᴰ*Lₗᵢₘ₁ᴰˢⁱ*min(5.4, ((4.4*exp(-4.23*Fₗᵢₘ₁ᴰˢⁱ)*Fₗᵢₘ₂ᴰˢⁱ + 1)*(1 + 2*Lₗᵢₘ₂ᴰˢⁱ)))   #22
end

#Diatoms have variable quadratic mortality term.
@inline function D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    return wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ) #eq13
end

#Phytoplankton forcing
@inline function (bgc::PISCES)(::Val{:P}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    δᴾ = bgc.exudation_of_DOC.P
    mᴾ = bgc.phytoplankton_mortality_rate.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    αᴾ = bgc.initial_slope_of_PI_curve.P

    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)

    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    #Grazing
    gₚᶻ = grazing_Z(P, D, POC, T, bgc)[2]     
    gₚᴹ =  grazing_M(P, D, Z, POC, T, bgc)[2]
    
    #Phytoplankton growth
    Lₗᵢₘᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    PARᴾ = P_PAR(PAR₁, PAR₂, PAR₃, bgc)
    μᴾ = phytoplankton_growth_rate(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)

    return (1-δᴾ)*μᴾ*P - mᴾ*concentration_limitation(P, Kₘ)*P - sh*wᴾ*P^2 - gₚᶻ*Z - gₚᴹ*M #eq 1
end

#Diatom forcing
@inline function (bgc::PISCES)(::Val{:D}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    δᴰ = bgc.exudation_of_DOC.D
    mᴰ = bgc.phytoplankton_mortality_rate.D
    Kₘ = bgc.half_saturation_const_for_mortality
    αᴰ = bgc.initial_slope_of_PI_curve.D

    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)

    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    #Grazing
    g_Dᶻ = grazing_Z(P, D, POC, T, bgc)[3]
    g_Dᴹ = grazing_M(P, D, Z, POC, T, bgc)[3]

    #Nutrient limitation
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]

    #Also required
    PARᴰ = D_PAR(PAR₁, PAR₂, PAR₃, bgc)
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)
    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc) #13

    return (1-δᴰ)*μᴰ*D - mᴰ*concentration_limitation(D, Kₘ)*D - sh*wᴰ*D^2 - g_Dᶻ*Z - g_Dᴹ*M    #eq 9
end

#Forcing for chlorophyll biomass of nanophytoplankton
@inline function (bgc::PISCES)(::Val{:Pᶜʰˡ}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
   #Parameters
    δᴾ = bgc.exudation_of_DOC.P
    αᴾ = bgc.initial_slope_of_PI_curve.P
    θₘᵢₙᶜʰˡ = bgc.min_ChlC_ratios_of_phytoplankton
    θₘₐₓᶜʰˡᴾ = bgc.max_ChlC_ratios_of_phytoplankton.P
    mᴾ = bgc.phytoplankton_mortality_rate.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton

    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)

    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    #Grazing
    gₚᶻ = grazing_Z(P, D, POC, T, bgc)[2]    
    gₚᴹ = grazing_M(P, D, Z, POC, T, bgc)[2]

    #Phytoplankton growth
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    Lₗᵢₘᴾ= P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    PARᴾ = P_PAR(PAR₁, PAR₂, PAR₃, bgc)
    μᴾ = phytoplankton_growth_rate(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)

    μ̌ᴾ = μᴾ / day_dependent_growth_rate(L_day) #15b
    ρᴾᶜʰˡ = 144*μ̌ᴾ * P / (αᴾ* Pᶜʰˡ* ((PARᴾ)/(L_day + eps(0.0))) + eps(0.0)) #15a

    return ((1-δᴾ)*(12*θₘᵢₙᶜʰˡ + (θₘₐₓᶜʰˡᴾ - θₘᵢₙᶜʰˡ)*ρᴾᶜʰˡ)*μᴾ*P - mᴾ*concentration_limitation(P, Kₘ)*Pᶜʰˡ 
             - sh*wᴾ*P*Pᶜʰˡ - nutrient_quota(Pᶜʰˡ, P)*gₚᶻ*Z - nutrient_quota(Pᶜʰˡ, P)*gₚᴹ*M)  #14
end

#Forcing for chlorophyll biomass of diatoms
@inline function (bgc::PISCES)(::Val{:Dᶜʰˡ}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    δᴰ = bgc.exudation_of_DOC.D
    αᴰ = bgc.initial_slope_of_PI_curve.D
    θₘᵢₙᶜʰˡ = bgc.min_ChlC_ratios_of_phytoplankton
    θₘₐₓᶜʰˡᴰ = bgc.max_ChlC_ratios_of_phytoplankton.D
    mᴰ = bgc.phytoplankton_mortality_rate.D
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms

    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)
    
    #Grazing
    g_Dᶻ = grazing_Z(P, D, POC, T, bgc)[3]
    g_Dᴹ = grazing_M(P, D, Z, POC, T, bgc)[3]

    #Diatom growth
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    PARᴰ = D_PAR(PAR₁, PAR₂, PAR₃, bgc)
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    μ̌ᴰ = μᴰ / (day_dependent_growth_rate(L_day) + eps(0.0)) #15b
    ρᴰᶜʰˡ = 144*μ̌ᴰ * D / (αᴰ* Dᶜʰˡ* ((PARᴰ)/(L_day + eps(0.0))) + eps(0.0)) #15a

    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc) #13
  
    return ((1-δᴰ)*(12*θₘᵢₙᶜʰˡ + (θₘₐₓᶜʰˡᴰ - θₘᵢₙᶜʰˡ)*ρᴰᶜʰˡ)*μᴰ*D
            - mᴰ*concentration_limitation(D, Kₘ)*Dᶜʰˡ - sh*wᴰ*D*Dᶜʰˡ 
            - nutrient_quota(Dᶜʰˡ, D)*g_Dᶻ*Z - nutrient_quota(Dᶜʰˡ, D)*g_Dᴹ*M)    #14
end

#Forcing for iron biomass of nanophytoplankton
@inline function (bgc::PISCES)(::Val{:Pᶠᵉ}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    δᴾ = bgc.exudation_of_DOC.P
    θₘₐₓᶠᵉᵖ = bgc.max_iron_quota.P
    mᴾ = bgc.phytoplankton_mortality_rate.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton.P
    K_Feᴾᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.P # this seems wrong as doesn't quite match parameter list
    Pₘₐₓ = bgc.threshold_concentration_for_size_dependency.P

    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)

    #Grazing
    gₚᶻ = grazing_Z(P, D, POC, T, bgc)[2]    
    gₚᴹ = grazing_M(P, D, Z, POC, T, bgc)[2]

    #Phytoplankton iron growth
    L_Feᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[6]
    bFe =  Fe   #defined in previous PISCES model
    μᴾᶠᵉ = phytoplankton_iron_biomass_growth_rate(P, Pᶠᵉ, θₘₐₓᶠᵉᵖ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe, T, bgc)

    return ((1-δᴾ)*μᴾᶠᵉ*P - mᴾ*concentration_limitation(P, Kₘ)*Pᶠᵉ - sh*wᴾ*P*Pᶠᵉ
             - nutrient_quota(Pᶠᵉ, P)*gₚᶻ*Z - nutrient_quota(Pᶠᵉ, P)*gₚᴹ*M ) #16
end

#Forcing for chlorophyll biomass of diatoms
@inline function (bgc::PISCES)(::Val{:Dᶠᵉ}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    δᴰ = bgc.exudation_of_DOC.D
    θₘₐₓᶠᵉᴰ = bgc.max_iron_quota.D
    mᴰ = bgc.phytoplankton_mortality_rate.D
    Kₘ = bgc.half_saturation_const_for_mortality
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton.D
    Dₘₐₓ = bgc.threshold_concentration_for_size_dependency.D
    K_Feᴰᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.D

    #Also required
    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc) #13
    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)
    
    #Limiting nutrients
    L = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    Lₗᵢₘᴰ = L[1]
    L_Feᴰ = L[6]

    #Grazing
    g_Dᶻ = grazing_Z(P, D, POC, T, bgc)[3]
    g_Dᴹ = grazing_M(P, D, Z, POC, T, bgc)[3]
   
    #Diatom iron growth
    bFe = Fe
    μᴰᶠᵉ = phytoplankton_iron_biomass_growth_rate(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Dₘₐₓ, L_Feᴰ, bFe, T, bgc)

    return ((1-δᴰ)*μᴰᶠᵉ*D - mᴰ*concentration_limitation(D, Kₘ)*Dᶠᵉ - sh*wᴰ*D*Dᶠᵉ
            - nutrient_quota(Dᶠᵉ, D)*g_Dᶻ*Z - nutrient_quota(Dᶠᵉ, D)*g_Dᴹ*M)    #16
end

#Forcing equations for silicon biomass of diatoms
@inline function (bgc::PISCES)(::Val{:Dˢⁱ}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)    #ϕ is latitude
    #Parameters
    δᴰ = bgc.exudation_of_DOC.D
    mᴰ = bgc.phytoplankton_mortality_rate.D
    Kₘ = bgc.half_saturation_const_for_mortality
    αᴰ = bgc.initial_slope_of_PI_curve.D

    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    #Grazing
    g_Dᶻ = grazing_Z(P, D, POC, T, bgc)[3]
    g_Dᴹ = grazing_M(P, D, Z, POC, T, bgc)[3]
 
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D

    #Diatom growth
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    PARᴰ = D_PAR(PAR₁, PAR₂, PAR₃, bgc)
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    #Also required
    θₒₚₜˢⁱᴰ = variation_in_SiC_ratio(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ, Si̅, bgc)
    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc) #13
    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)
    
    return (θₒₚₜˢⁱᴰ*(1-δᴰ)*μᴰ*D - nutrient_quota(Dˢⁱ, D)*g_Dᴹ*M -  nutrient_quota(Dˢⁱ, D)*g_Dᶻ*Z
           - mᴰ*concentration_limitation(D, Kₘ)*Dˢⁱ - sh*wᴰ*D*Dˢⁱ) #21
end