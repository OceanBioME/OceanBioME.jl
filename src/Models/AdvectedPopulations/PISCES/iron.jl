@inline function free_iron(::SimpleIron, Fe, DOC, T)
    # maybe some of these numbers should be parameters
    ligands = max(0.6, 0.09 * (DOC + 40) - 3)
    K = exp(16.27 - 1565.7 / max(T + 273.15, 5))
    Δ = 1 + K * ligands - K * Fe

    return (-Δ + √(Δ^2 + 4K * Fe)) / 2K
end

# this should be dispatched on an abstract type if we implement complex chemistry
@inline function aggregation_of_colloidal_iron(iron::SimpleIron, dom, bgc, z, DOC, POC, GOC, Fe, T, zₘₓₗ)
    _, Φ₁, Φ₂, Φ₃ = aggregation(dom, bgc, z, DOC, POC, GOC, zₘₓₗ)

    Fe′ = free_iron(iron, Fe, DOC, T)
    ligand_iron = Fe - Fe′
    colloidal_iron = 0.5 * ligand_iron

    CgFe1 = (Φ₁ + Φ₃) * colloidal_iron / DOC
    CgFe2 = Φ₂ * colloidal_iron / DOC

    return CgFe1 + CgFe2, CgFe1, CgFe2
end

# This document contains functions for the following:
    # free_iron(eq65), dissolved free inorganic iron
    # iron_colloid_aggregation_1, iron_colloid_aggregation_2, enhanced_scavenging, bacterial_uptake_Fe (eqs 61, 62, 63)
    # Forcing for Fe (eq60)


#Iron is modelled in PISCES using simple chemistry model. Iron is assumed to be in the form of free organic iron Fe', and complexed iron FeL.
#Iron is explictly modelled in the following compartments of the model, Pᶠᵉ, Dᶠᵉ, SFe, BFe, Fe. In zooplankton a fixed iron carbon ratio is assumed, and iron implictly modelled in zooplankton in this way. Iron in bacteria is formulated through Bactfe.
#Iron is lost through Scav and enhanced_scavenging terms. Aside from this, iron is conserved, accounting for all other explicit and implicit compartments of iron

#Determine concentration of free organic iron. This is the only form of iron assumed susceptible to scavenging.
@inline function free_iron(Fe, DOC, T)
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6) # This may also be taken to be a constant parameter, total_concentration_of_iron_ligands
    K_eqᶠᵉ = exp(16.27 - 1565.7/max(T + 273.15, 5)) #check this value
    Δ = 1 +  K_eqᶠᵉ*Lₜ -  K_eqᶠᵉ*Fe
    return (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ*Fe))/(2*K_eqᶠᵉ + eps(0.0)) #eq65
end

#Colloids of iron may aggregate with DOM and be transferred to particulate pool
#iron_colloid_aggregation_1 is aggregation of colloids with DOC and POC. Routed to SFe.
@inline function iron_colloid_aggregation_1(sh, Fe, POC, DOC, T, bgc)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
   
    FeL = Fe - free_iron(Fe, DOC, T) #eq64
    Fe_coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll #eq61a
end

#iron_colloid_aggregation_2 is aggregation of colloids with GOC. Routed to BFe.
@inline function iron_colloid_aggregation_2(sh, Fe, T, DOC, GOC, bgc)
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    FeL = Fe - free_iron(Fe, DOC, T)
    Fe_coll = 0.5*FeL
    return a₃*GOC*sh*Fe_coll #eq61b
end

#When dissolved iron concentrations exceed total ligand concentrations scavenging is enhanced.
@inline function enhanced_scavenging(Fe, DOC, T, bgc)
    λᶠᵉ = 1e-3 * bgc.slope_of_scavenging_rate_of_iron #parameter not defined in parameter list. Assumed scaled version λ_Fe to fit dimensions of Fe¹.
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6)
    return λᶠᵉ*max(0, Fe - Lₜ)*free_iron(Fe, DOC, T) #eq62
end

#Formulation for bacterial uptake of iron.
@inline function bacterial_uptake_Fe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)
    K_Feᴮ¹ = bgc.Fe_half_saturation_const_for_Bacteria
    θₘₐₓᶠᵉᵇᵃᶜᵗ = bgc.max_FeC_ratio_of_bacteria
    Bact = bacterial_biomass(zₘₐₓ, z, Z, M) 
    Lₗᵢₘᵇᵃᶜᵗ = bacterial_activity(DOC, PO₄, NO₃, NH₄, bFe, bgc)[2]
    bₚ = bgc.temperature_sensitivity_of_growth
    return μₘₐₓ⁰*(bₚ^T)*Lₗᵢₘᵇᵃᶜᵗ*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe + eps(0.0)) #eq63
end

@inline function (bgc::PISCES)(::Val{:Fe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃) #eq60
    #Parameters
    σᶻ, σᴹ = bgc.non_assimilated_fraction
    eₘₐₓᶻ, eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton
    δᴾ, δᴰ = bgc.exudation_of_DOC
    θₘₐₓᶠᵉᴾ, θₘₐₓᶠᵉᴰ = bgc.max_iron_quota
    Sᵣₐₜᴾ, Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton
    K_Feᴾᶠᵉᵐⁱⁿ, K_Feᴰᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake
    Pₘₐₓ, Dₘₐₓ = bgc.threshold_concentration_for_size_dependency
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M
    wₚₒ = bgc.sinking_speed_of_POC

    γᴹ = bgc.excretion_as_DOM.M #Removed γᴹ factor from upper_respiration to conserve iron implicitly lost through mesozooplankton quadratic mortality.

    bFe = Fe

    #Growth rate of iron biomass of phytoplankton
    L_Feᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[6]
    L_Feᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[6]

    μᴾᶠᵉ = phytoplankton_iron_biomass_growth_rate(P, Pᶠᵉ, θₘₐₓᶠᵉᴾ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe, T, bgc)
    μᴰᶠᵉ = phytoplankton_iron_biomass_growth_rate(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Dₘₐₓ, L_Feᴰ, bFe, T, bgc)

    #Iron quotas
    θᶠᵉᴾ = nutrient_quota(Pᶠᵉ, P)
    θᶠᵉᴰ = nutrient_quota(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = nutrient_quota(SFe, POC)
    θᶠᵉᴳᴼᶜ = nutrient_quota(BFe, GOC)
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = grazing_M(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ, gₚₒ_FF, g_GOC_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    w_GOC = sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*gₚᶻ + θᶠᵉᴰ*g_Dᶻ + θᶠᵉᴾᴼᶜ*gₚₒᶻ #over P, D, POC
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*gₚᴹ + θᶠᵉᴰ*g_Dᴹ + θᶠᵉᴾᴼᶜ*gₚₒᴹ + θᶠᵉᶻ*g_Zᴹ #graze on P, D, POC, Z 

    #Iron in bacteria
    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ))
    Bactfe = bacterial_uptake_Fe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)

    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc) #eₘₐₓᶻ used in paper but changed here to be consistent with eqs 24, 28
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)
    λₚₒ¹ = particles_carbon_degradation_rate(T, O₂, bgc)

    return (max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/(∑gᶻ + eps(0.0))) - eᶻ*θᶠᵉᶻ)*∑gᶻ*Z 
            + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ )/(∑gᴹ+∑g_FFᴹ + eps(0.0)) - eᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M 
            + θᶠᵉᶻ*upper_respiration(M, T, bgc) + λₚₒ¹*SFe 
            - (1 - δᴾ)*μᴾᶠᵉ*P - (1 - δᴰ)*μᴰᶠᵉ*D 
            - Fe_scavenging(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) - iron_colloid_aggregation_1(sh, Fe, POC, DOC, T, bgc)
            - iron_colloid_aggregation_2(sh, Fe, T, DOC, GOC, bgc) - enhanced_scavenging(Fe, DOC, T, bgc) - Bactfe) #eq60
end
