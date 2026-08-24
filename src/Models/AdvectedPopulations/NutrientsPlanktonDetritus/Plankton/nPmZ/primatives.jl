@inline temperature_scaling(Q10, T₀, T) = Q10 ^ ((T - T₀) / oftype(T, 10))

@inline eppley_temperature_scaling(T) =
    oftype(T, 0.12) * max(zero(T), min(T, oftype(T, 27))) ^ oftype(T, 0.4)

@inline autotroph_temperature_scaling(a::Autotroph{traits}, T) where traits =
    _autotroph_temperature_scaling(Val(traits.temperature_function), a, T)

@inline _autotroph_temperature_scaling(::Val{:q10}, a, T) =
    temperature_scaling(a.temperature_sensitivity, a.temperature_reference, T)

@inline _autotroph_temperature_scaling(::Val{:power}, a, T) = eppley_temperature_scaling(T)

@inline max_specific_rate(PCref, nutrient_limitation, Tf, T, temp_thres) = (PCref * nutrient_limitation * Tf) * (T >= temp_thres)

@subcolumn_average @inline function light_limited_rate(PCmax, α, θᶜ, PAR, ε)
    x = α * θᶜ * PAR
    return PCmax * (1 - exp(-x / (PCmax + ε)))
end


# Liebig nutrient_limitation over the nutrients that apply to a PFT (silicifier/nfixer baked at gen time).
# DOP = 0 in Phase 3 ⇒ VPtot = VPO₄. VSi divides by (SiO₃+kSi); for non-silicifiers it is selected
# out by `ifelse` (never ×Bool — the dead branch can be NaN when kSi = 0 and SiO₃ = 0).
@inline competitive_limitation(a, b, ka, kb) = (a/ka + b/kb) / (1 + a/ka + b/kb)

@inline function minimum_limitation(silicifier, nfixer, kNO₃, kNH₄, kPO₄, kDOP, kFe, kSi,
                                    NO₃, NH₄, PO₄, DOP, Fe, SiO₃)
    # isn't there evidence that sometimes nitrifiers also do uptake
    VN  = ifelse(nfixer, one(NO₃), competitive_limitation(NO₃, max(zero(NH₄), NH₄), kNO₃, kNH₄))
    VFe = Fe / (Fe + kFe)
    VP  = competitive_limitation(PO₄, DOP, kPO₄, kDOP)

    f   = min(VN, VFe, VP)

    VSi = SiO₃ / (SiO₃ + kSi)

    return ifelse(silicifier, min(f, VSi), f)
end

# PO₄-only Monod. The picpoc P-limitation term uses THIS, not the competitive VP total.
@inline phosphate_limitation(PO₄, DOP, kPO₄, kDOP) =
    (PO₄ / kPO₄) / (1 + PO₄ / kPO₄ + DOP / kDOP)

@inline carbon_limitation(::Val{false}, f, a, p, i, j, k) = f
@inline function carbon_limitation(::Val{true}, f, a, p, i, j, k)
    @inbounds CO₂ = p.carbon_dioxide[i, j, k]
    VCO₂ = CO₂ / (CO₂ + a.carbon_dioxide_half_saturation)
    return min(f, VCO₂)
end

# Galbraith & Martiny growth P:C, capped at min N:P. MARBL gQp.
@inline growth_phosphorus_quota(sₚ, iₚ, Qp_max, PO₄::FT) where FT = 
    min((sₚ * PO₄ + iₚ) * convert(FT, 1e-3), Qp_max)

# MARBL gQfe.
@inline function growth_iron_quota(gQfe₀, gQfe_min, r_Fe, kFe, Fe)
    threshold = r_Fe * kFe
    economised = max(gQfe₀ * Fe / threshold, gQfe_min)
    return ifelse(Fe < threshold, economised, gQfe₀)
end

@inline function growth_silicon_quota(gQsi₀, gQsi_max, gQsi_min, r_Fe, r_Si, kFe, kSi, Fe, SiO₃)
    Fe_threshold = r_Fe * kFe
    Si_threshold = r_Si * kSi

    low_Fe = (Fe > 0) & (Fe < Fe_threshold) & (SiO₃ > Si_threshold)
    raised = ifelse(Fe == 0, gQsi_max,
                    ifelse(low_Fe, min(gQsi₀ * Fe_threshold / Fe, gQsi_max), gQsi₀))

    lowered = max(raised * SiO₃ / Si_threshold, gQsi_min)
    return ifelse(SiO₃ < Si_threshold, lowered, raised)
end

# GD98 dynamic Chl:C
@subcolumn_average @inline function photoacclimation_rate(a, p, PCmax, θᶜ, Chl, PAR)
    θᴺmax = a.maximum_chlorophyll_to_nitrogen
    α = a.initial_PI_slope
    Q = p.nitrogen_to_carbon

    PCphoto = light_limited_rate(PCmax, a.initial_PI_slope, θᶜ, PAR, p.growth_rate_regularisation)
    x = α * θᶜ * PAR
    ρᶜʰˡ = θᴺmax * PCphoto / x
    photoacc = ρᶜʰˡ * PCphoto * Q / θᶜ * Chl

    ifelse((x > 0) & (θᶜ > 0), photoacc, zero(photoacc))
end

@inline function aggregation_rate(mort2, mort2_exp, agg_min, agg_max, P′)
    agg = min(agg_max * P′, mort2 * P′ ^ mort2_exp)
    return max(agg_min * P′, agg)
end

@inline function calcite_formation(f_prod, f_photo, C_thres, Tc1, Tc2, photosynthesis, nutrient_limitation, T, C)
    cf = f_prod * photosynthesis * nutrient_limitation * nutrient_limitation
    cf = ifelse(T < Tc1, cf * max(T - Tc2, zero(T)) / (Tc1 - Tc2), cf)
    return ifelse(C > C_thres, min(cf * C / C_thres, f_photo * photosynthesis), cf)
end

# for coccolithophores Krumhardt et al. 2019/2020
@inline function picpoc_calcite_formation(CO₂, VPO₄, photosynthesis, T)
    ϖ = min(one(T), max(zero(T), oftype(T, 1e-4) * T^3))            # temperature effect (cubic)
    ϖ = max(zero(T), oftype(T, -0.0136) * CO₂ + ϖ + oftype(T, 0.21)) # CO₂ effect (picpoc ↓ as CO₂ ↑)
    ϖ = max(zero(T), oftype(T, -0.48) * VPO₄ + ϖ + oftype(T, 0.48))  # P-limitation (picpoc ↑ as VPO₄ ↓)
    return ϖ * photosynthesis
end

@inline nitrogen_fixation(r_Nfix, Q, photosynthesis, NO₃_V, NH₄_V) = Q * photosynthesis * r_Nfix - NO₃_V - NH₄_V
@inline nitrogen_excretion(r_Nfix, Q, photosynthesis, NO₃_V, NH₄_V) =
    nitrogen_fixation(r_Nfix, Q, photosynthesis, NO₃_V, NH₄_V) + NO₃_V + NH₄_V - Q * photosynthesis

# grazing rate on a prey class of biomass B by a predator of biomass Cz
@inline function grazing_rate(umax, kgrz, sigmoidal, Tf_z, Cz, B)
    saturation = ifelse(sigmoidal, B * B / (B * B + kgrz * kgrz), B / (B + kgrz))
    return umax * Tf_z * Cz * saturation * (B > 0)
end

# ballast enhancement of grazed-→-POC routing
@inline calcifier_graze_poc_fraction(CaCO₃_poc_min, QCaCO₃, spc_poc_fac, P′, f_lim) =
    max(CaCO₃_poc_min * QCaCO₃,
        min(spc_poc_fac * (P′ + oftype(P′, 0.6)) ^ oftype(P′, 1.6), f_lim))

@inline zoo_loss_rate(mort2, mort, mort2_exp, Z′, Tf_z) = (mort2 * Z′ ^ mort2_exp + mort * Z′) * Tf_z

