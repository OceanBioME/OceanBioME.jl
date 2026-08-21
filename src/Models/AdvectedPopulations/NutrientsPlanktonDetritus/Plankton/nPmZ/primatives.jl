@inline temperature_scaling(Q10, T₀, T) = Q10 ^ ((T - T₀) / 10)

@inline eppley_temperature_scaling(T) = 0.12 * max(zero(T), min(T, oftype(T, 27))) ^ oftype(T, 0.4)

@inline autotroph_temperature_scaling(a::Autotroph{traits}, T) where traits =
    _autotroph_temperature_scaling(Val(traits.temperature_function), a, T)

@inline _autotroph_temperature_scaling(::Val{:q10}, a, T) =
    temperature_scaling(a.temperature_sensitivity, a.temperature_reference, T)

@inline _autotroph_temperature_scaling(::Val{:power}, a, T) = eppley_temperature_scaling(T)

@inline max_specific_rate(PCref, nutrient_limitation, Tf, T, temp_thres) = (PCref * nutrient_limitation * Tf) * (T >= temp_thres)

@inline function light_limited_rate(PCmax, α, θᶜ, PAR, ε)
    x = α * θᶜ * PAR
    return PCmax * (1 - exp(-x / (PCmax + ε)))
end

# PAR sub column _graze_weighted_detritus
# the return value is the area fraction φⱼ and shortwave ratio rⱼ
# where PARⱼ(z) = PAR(k)⋅rⱼ, but every non-linear function, is Σⱼ φⱼ·f(PAR·rⱼ)
# This is to allow sub-column shading from ice (but this isn't implemented in ClimaSeaIce yet)
# when no sub-columns are specified return (1, 1)
@inline par_subcolumns(i, j, aux::NamedTuple{names}) where {names} = _par_subcolumns(i, j, aux, Val(:PAR_fractions in names))

@inline _par_subcolumns(i, j, aux, ::Val{true})  = @inbounds (aux.PAR_areas[i, j, 1], aux.PAR_fractions[i, j, 1])
@inline _par_subcolumns(i, j, aux, ::Val{false}) = (one(defaults.FloatType), one(defaults.FloatType))

# sometimes we need the light at the interface between two cells rather than the cell
# center so we can either pass it as an auxiliary field (for the testing) or interpolate
@inline interface_par(i, j, k, grid, aux::NamedTuple{names}) where {names} =
    _interface_par(i, j, k, grid, aux, Val(:PAR_interface in names))

@inline _interface_par(i, j, k, grid, aux, ::Val{true})  = @inbounds aux.PAR_interface[i, j, k]
@inline _interface_par(i, j, k, grid, aux, ::Val{false}) = ℑzᵃᵃᶠ(i, j, k, grid, aux.PAR)

# Liebig nutrient_limitation over the nutrients that apply to a PFT (silicifier/nfixer baked at gen time).
# DOP = 0 in Phase 3 ⇒ VPtot = VPO₄. VSi divides by (SiO₃+kSi); for non-silicifiers it is selected
# out by `ifelse` (never ×Bool — the dead branch can be NaN when kSi = 0 and SiO₃ = 0).
@inline competitive_limitation(a, b, ka, kb) = (a/ka + b/kb) / (1 + a/ka + b/kb)

@inline function minimum_limitation(silicifier, nfixer, kNO₃, kNH₄, kPO₄, kDOP, kFe, kSi,
                                    NO₃, NH₄, PO₄, DOP, Fe, SiO₃)
    # isn't there evidence that sometimes nitrifiers also do uptake
    VN  = ifelse(nfixer, one(a), competitive_limitation(NO₃, max(zero(NH₄), NH₄), kNO₃, kNH₄))
    VFe = Fe / (Fe + kFe)
    VP  = competitive_limitation(PO₄, DOP, kPO₄, kDOP)

    f   = min(VN, VFe, VP)

    VSi = SiO₃ / (SiO₃ + kSi)

    return ifelse(silicifier, min(f, VSi), f)
end

@inline carbon_limitation(::Val{false}, f, a, p, i, j, k) = f
@inline function carbon_limitation(::Val{true}, f, a, p, i, j, k)
    @inbounds CO₂ = p.carbon_dioxide[i, j, k]
    VCO2 = CO₂ / (CO₂ + a.carbon_dioxide_half_saturation)
    return min(f, VCO2)
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
@inline function photoacclimation_rate(PCphoto, θᴺmax, Q, θᶜ, Chl, x)
    ρᶜʰˡ = θᴺmax * PCphoto / x
    photoacc = ρᶜʰˡ * PCphoto * Q / θᶜ * Chl
    return ifelse((x > 0) & (θᶜ > 0), photoacc, zero(photoacc))
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
    π = min(one(T), max(zero(T), 1e-4 * T^3))          # temperature effect (cubic)
    π = max(zero(T), -0.0136 * CO₂  + π + 0.21)          # CO₂ effect (picpoc ↓ as CO₂ ↑)
    π = max(zero(T), -0.48  * VPO₄ + π + 0.48)           # P-limitation effect (picpoc ↑ as VPO₄ ↓)
    return π * photosynthesis
end

@inline nitrogen_fixation(r_Nfix, Q, photosynthesis, NO₃_V, NH4_V) = Q * photosynthesis * r_Nfix - NO₃_V - NH4_V
@inline nitrogen_excretion(r_Nfix, Q, photosynthesis, NO₃_V, NH4_V) =
    nitrogen_fixation(r_Nfix, Q, photosynthesis, NO₃_V, NH4_V) + NO₃_V + NH4_V - Q * photosynthesis

# grazing rate on a prey class of biomass B by a predator of biomass Cz
@inline function grazing_rate(umax, kgrz, sigmoidal, Tf_z, Cz, B)
    saturation = ifelse(sigmoidal, B * B / (B * B + kgrz * kgrz), B / (B + kgrz))
    return umax * Tf_z * Cz * saturation * (B > 0)
end

# ballast enhancement of grazed-→-POC routing
@inline calcifier_graze_poc_fraction(CaCO₃_poc_min, QCaCO₃, spc_poc_fac, P′, f_lim) =
    max(CaCO₃_poc_min * QCaCO₃, min(spc_poc_fac * (P′ + 0.6) ^ 1.6, f_lim))

@inline zoo_loss_rate(mort2, mort, mort2_exp, Z′, Tf_z) = (mort2 * Z′ ^ mort2_exp + mort * Z′) * Tf_z

