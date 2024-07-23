
#This document contains functions for:
    #μₙₒ₃ᴶ,  μₙₕ₄ᴶ (eq8)
    #ΔO₂ (eq57)
    #Nitrif (eq56)
    #N_fix (eq58)
    #Forcing for NO₃ and NH₄ (eqs54, 55)

@inline function get_μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T,  zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc) 
    αᴾ = bgc.initial_slope_of_PI_curve.P
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[3]
    return μᴾ * K_mondo(Lₙₒ₃ᴾ, Lₙₕ₄ᴾ) #eq8
end

@inline function get_μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    αᴾ = bgc.initial_slope_of_PI_curve.P
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[3]
    return μᴾ * K_mondo(Lₙₕ₄ᴾ, Lₙₒ₃ᴾ) #eq8
end

@inline function get_μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc) 
    αᴰ = bgc.initial_slope_of_PI_curve.D
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴰ =  μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)
    Lₙₒ₃ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[4]
    Lₙₕ₄ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[3]
    return μᴰ * K_mondo(Lₙₒ₃ᴰ, Lₙₕ₄ᴰ) #eq8
end

@inline function get_μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    αᴰ = bgc.initial_slope_of_PI_curve.D
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴰ =  μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)
    Lₙₒ₃ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[4]
    Lₙₕ₄ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[3]
    return μᴰ * K_mondo(Lₙₕ₄ᴰ, Lₙₒ₃ᴰ) #eq8
end

@inline function ΔO₂(O₂, bgc)
    O₂ᵐⁱⁿ¹ = bgc.half_sat_const_for_denitrification1
    O₂ᵐⁱⁿ² = bgc.half_sat_const_for_denitrification2

    return min(1, max(0.4*(O₂ᵐⁱⁿ¹-O₂)/(O₂ᵐⁱⁿ²+O₂))) #eq57
end

@inline Nitrif(NH₄, O₂, λₙₕ₄, PAR, bgc) = λₙₕ₄*NH₄*(1-ΔO₂(O₂, bgc))/(1+PAR) #eq56a

@inline function (bgc::PISCES)(::Val{:NO₃}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 

    λₙₕ₄ =  bgc.max_nitrification_rate

    Rₙₕ₄ = bgc.NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER
    Rₙₒ₃ = bgc.NC_stoichiometric_ratio_of_dentitrification

    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)

    bFe = Fe
    zₘₐₓ = max(zₑᵤ, zₘₓₗ) #35a
    Bact = get_Bact(zₘₐₓ, z, Z, M)

    μₙₒ₃ᴾ = get_μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₒ₃ᴰ = get_μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)

    return Nitrif(NH₄, O₂, λₙₕ₄, PAR, bgc) - μₙₒ₃ᴾ*P - μₙₒ₃ᴰ*D - Rₙₕ₄*λₙₕ₄*ΔO₂(O₂, bgc)*NH₄ - Rₙₒ₃*get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc)
end

# The following relate specifically to NH₄ forcing
#Change to ifelse

@inline Lₙᴰᶻ(Lₙᴾ) = ifelse(Lₙᴾ>=0.08, 0.01, 1 - Lₙᴾ) #eq58
    
@inline function N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc) #eq 58b
    N_fixᵐ = bgc.max_rate_of_nitrogen_fixation
    K_Feᴰᶻ = bgc.Fe_half_saturation_constant_of_nitrogen_fixation
    Kₚₒ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate.P
    E_fix = bgc.photosynthetic_parameter_of_nitrogen_fixation
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    μₚ = μ⁰ₘₐₓ*fₚ(T, bgc)
    Lₙᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[5]

    return N_fixᵐ*max(0,μₚ - 2.15)*Lₙᴰᶻ(Lₙᴾ)*min(K_mondo(bFe, K_Feᴰᶻ), K_mondo(PO₄, Kₚₒ₄ᴾᵐⁱⁿ))*(1 - exp((-PAR/E_fix)))
end

@inline function (bgc::PISCES)(::Val{:NH₄}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 
    
    #Parameters
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M

    bFe = Fe #Change this!

    #L_day
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc) 
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = get_grazingᴹ(P, D, Z, POC, T, bgc) 

    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
        #g_Z error
    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    μₙₕ₄ᴾ = get_μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₕ₄ᴰ = get_μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)

    zₘₐₓ = max(zₑᵤ, zₘₓₗ) #35a
    Bact = get_Bact(zₘₐₓ, z, Z, M)
   
    return γᶻ*(1-eᶻ-σᶻ)*∑gᶻ*Z + γᴹ*(1-eᴹ-σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*Rᵤₚ(M, T, bgc) + get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) + N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc) - Nitrif(NH₄, O₂, λₙₕ₄, PAR, bgc) - λₙₕ₄*ΔO₂(O₂, bgc)*NH₄ - μₙₕ₄ᴾ*P - μₙₕ₄ᴰ*D
end
