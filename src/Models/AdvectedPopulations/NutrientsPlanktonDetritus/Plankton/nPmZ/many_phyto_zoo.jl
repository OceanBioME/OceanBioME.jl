#####
##### Combined class
#####
struct ManyPhytoZoo{LN, P, Z, AT, ZT, CD, FT} <: AbstractPlankton{LN}
                                  autotrophs :: AT  # NamedTuple of Autotroph
                                 zooplankton :: ZT  # NamedTuple of Heterotroph
                              carbon_dioxide :: CD  # mmol/m³; aqueous CO₂ (nothing if unused)

             zooplankton_phosphate_to_carbon :: FT  # mmol P/mmol C; MARBL Qp_zoo
                  zooplankton_iron_to_carbon :: FT  # mmol Fe/mmol C; MARBL Qfe_zoo

                          nitrogen_to_carbon :: FT  # mmol N/mmol C; MARBL Q
                  growth_rate_regularisation :: FT  # 1/s; MARBL epsTinv 
                concentration_regularisation :: FT  # mmol C/m³; MARBL epsC 
                             labile_fraction :: FT  # ; MARBL parm_labile_ratio
                    nitrogen_to_DON_fraction :: FT  # ; MARBL f_toDON (DON:DOC production scaling)
                  phosphorus_to_DOP_fraction :: FT  # ; MARBL f_toDOP (surplus autotroph P split to DOP vs PO₄)
                        aggregation_exponent :: FT  # ; MARBL auto_mort2_exp

                       phosphate_quota_slope :: FT  # ; MARBL PquotaSlope (Galbraith–Martiny; ×10⁻³ with [PO₄] in mmol m⁻³)
                   phosphate_quota_intercept :: FT  # ; MARBL PquotaIntercept (Galbraith–Martiny; ×10⁻³)
                 phosphate_to_carbon_maximum :: FT  # mol P/mol C; MARBL PquotaMinNP

                 iron_quota_threshold_factor :: FT  # ; MARBL gQ_Fe_kFe_thres (r_Fe)

                     silicon_quota_reference :: FT  # mol Si/mol C; MARBL gQsi_0
                       maximum_silicon_quota :: FT  # mol Si/mol C; MARBL gQsi_max
                       minimum_silicon_quota :: FT  # mol Si/mol C; MARBL gQsi_min
              silicon_quota_threshold_factor :: FT  # ; MARBL gQ_Si_kSi_thres (r_Si)

                 calcite_production_fraction :: FT  # ; MARBL parm_f_prod_sp_CaCO3
                 maximum_calcifying_fraction :: FT  # ; MARBL f_photosp_CaCO3
                     calcite_bloom_threshold :: FT  # mmol C/m³; MARBL CaCO3_sp_thres
               calcite_temperature_threshold :: FT  # °C; MARBL CaCO3_temp_thres1
                 calcite_temperature_minimum :: FT  # °C; MARBL CaCO3_temp_thres2
                       maximum_calcite_quota :: FT  # mol CaCO₃/mol C; MARBL QCaCO3_max

                     nitrogen_fixation_ratio :: FT  # ; MARBL r_Nfix_photo

                   loss_threshold_full_depth :: FT  # m; MARBL thres_z1_auto
                   loss_threshold_zero_depth :: FT  # m; MARBL thres_z2_auto

                    zoo_aggregation_exponent :: FT  # ; MARBL zoo_mort2_exp
               zoo_loss_threshold_full_depth :: FT  # m; MARBL thres_z1_zoo
               zoo_loss_threshold_zero_depth :: FT  # m; MARBL thres_z2_zoo

    grazed_calcite_remineralisation_fraction :: FT  # ; MARBL f_graze_CaCO3_remin
     grazed_silica_remineralisation_fraction :: FT  # ; MARBL f_graze_si_remin
                     calcite_ballast_minimum :: FT  # ; MARBL caco3_poc_min
                      small_phyto_poc_factor :: FT  # 1/(mmol C); MARBL spc_poc_fac
                grazed_small_phyto_poc_limit :: FT  # ; MARBL f_graze_sp_poc_lim

                  alkalinity_phosphate_ratio :: FT  # P:C used ONLY for the alkalinity PO₄ term (0 = MARBL §20.4, which
                                                    # omits it; set ≈1/117 to include the arguably-correct PO₄–Alk term)
end

const ManyPhytoZoo_NPD = NutrientsPlanktonDetritus{<:Any, <:Nutrients, <:ManyPhytoZoo, <:Any, <:AbstractInorganicCarbon}

const _manifested_marbl_plankton = Set{Tuple}()

@inline autotroph(p::ManyPhytoZoo, ::Val{S}) where S = getproperty(p.autotrophs, S)
@inline traits(p::ManyPhytoZoo, val_sname) = traits(autotroph(p, val_sname))


@inline chlorophyll_to_carbon_ratio(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    @inbounds autotroph_chlorophyll(val_prefix, i, j, k, p, fields) /
              (autotroph_carbon(val_prefix, i, j, k, p, fields) + p.concentration_regularisation)

@inline standing_phosphorus_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    @inbounds autotroph_phosphorus(val_prefix, i, j, k, p, fields) /
              (autotroph_carbon(val_prefix, i, j, k, p, fields) + p.concentration_regularisation)

@inline standing_iron_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    @inbounds autotroph_iron(val_prefix, i, j, k, p, fields) /
              (autotroph_carbon(val_prefix, i, j, k, p, fields) + p.concentration_regularisation)

@inline function standing_silicon_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    @inbounds q = autotroph_silicon(val_prefix, i, j, k, p, fields) /
                  (autotroph_carbon(val_prefix, i, j, k, p, fields) + p.concentration_regularisation)
    return min(q, p.maximum_silicon_quota)
end

@inline function standing_calcite_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    @inbounds q = autotroph_calcite(val_prefix, i, j, k, p, fields) /
                  (autotroph_carbon(val_prefix, i, j, k, p, fields) + p.concentration_regularisation)
    return min(q, p.maximum_calcite_quota)
end

@inline function growth_phosphorus_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    @inbounds PO₄ = fields.PO₄[i, j, k]

    return growth_phosphorus_quota(p.phosphate_quota_slope, p.phosphate_quota_intercept,
                                   p.phosphate_to_carbon_maximum, PO₄)
end

@inline function growth_iron_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)
    @inbounds Fe = fields.Fe[i, j, k]
    return growth_iron_quota(a.growth_iron_quota_reference, a.minimum_iron_quota,
                             p.iron_quota_threshold_factor, a.nutrient_half_saturations.iron, Fe)
end

@inline function growth_silicon_quota(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)
    @inbounds begin
        Fe = fields.Fe[i, j, k]; SiO₃ = fields.Si[i, j, k]
    end
    return growth_silicon_quota(p.silicon_quota_reference, p.maximum_silicon_quota,
                                p.minimum_silicon_quota, p.iron_quota_threshold_factor,
                                p.silicon_quota_threshold_factor,
                                a.nutrient_half_saturations.iron,
                                a.nutrient_half_saturations.silicate, Fe, SiO₃)
end

@inline function nutrient_limitation(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a   = autotroph(p, val_prefix)
    khs = a.nutrient_half_saturations

    @inbounds begin
        NO₃  = fields.NO₃[i, j, k]
        NH₄  = fields.NH₄[i, j, k]
        PO₄  = fields.PO₄[i, j, k]
        Fe   =  fields.Fe[i, j, k]
        DOP  = fields.DOP[i, j, k]
        SiO₃ =  fields.Si[i, j, k]
    end

    f = minimum_limitation(traits(a).silicifier, traits(a).nitrogen_fixer,
                           khs.nitrate, khs.ammonia, khs.phosphate, khs.DOP, khs.iron, khs.silicate,
                           NO₃, NH₄, PO₄, DOP, Fe, SiO₃)

    return carbon_limitation(Val(traits(a).carbon_limited), f, a, p, i, j, k)
end

#####
##### Photosynthesis and photoacclimation
#####

# gross carbon fixation (mmol C/m³/s), φ-weighted over the ice-radiation PAR sub-columns
@inline function photosynthesis(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields, aux)
    a = autotroph(p, val_prefix)

    @inbounds begin
        T   = fields.T[i, j, k]
        PAR = @preserve_subcolumns aux.PAR[i, j, k]
        C   = autotroph_carbon(val_prefix, i, j, k, p, fields)
    end

    fₙ    = nutrient_limitation(val_prefix, i, j, k, grid, p, fields)
    Tf    = autotroph_temperature_scaling(a, T)
    PCmax = max_specific_rate(a.photosynthesis_rate_reference, fₙ, Tf, T, a.temperature_threshold)
    θᶜ    = chlorophyll_to_carbon_ratio(val_prefix, i, j, k, grid, p, fields)

    PCphoto = light_limited_rate(PCmax, a.initial_PI_slope, θᶜ, PAR, p.growth_rate_regularisation)

    return PCphoto * C
end

# chlorophyll synthesis (mg Chl/m³/s), Geider et al. 1998 dynamic Chl:C
@inline function photoacclimation(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields, aux)
    a = autotroph(p, val_prefix)

    @inbounds begin
        T   = fields.T[i, j, k]
        PAR = @preserve_subcolumns aux.PAR[i, j, k]
        Chl = autotroph_chlorophyll(val_prefix, i, j, k, p, fields)
    end

    fₙ    = nutrient_limitation(val_prefix, i, j, k, grid, p, fields)
    Tf    = autotroph_temperature_scaling(a, T)
    PCmax = max_specific_rate(a.photosynthesis_rate_reference, fₙ, Tf, T, a.temperature_threshold)
    θᶜ    = chlorophyll_to_carbon_ratio(val_prefix, i, j, k, grid, p, fields)

    return photoacclimation_rate(a, p, PCmax, θᶜ, Chl, PAR)
end

#####
##### Autotroph losses
#####

# loss-active biomass P′, with the depth ramp that switches the loss threshold off at depth
@inline function loss_active_biomass(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)

    @inbounds begin
        T = fields.T[i, j, k]
        C = autotroph_carbon(val_prefix, i, j, k, p, fields)
    end

    depth = -znode(i, j, k, grid, Center(), Center(), Center())

    f_thres = clamp((p.loss_threshold_zero_depth - depth) /
                    (p.loss_threshold_zero_depth - p.loss_threshold_full_depth),
                    zero(depth), one(depth))

    thres = ifelse(T < a.temperature_threshold,
                   a.cold_loss_threshold_concentration, a.loss_threshold_concentration)

    return max(C - f_thres * thres, zero(C))
end

@inline function autotroph_linear_loss(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)

    @inbounds T = fields.T[i, j, k]

    P′ = loss_active_biomass(val_prefix, i, j, k, grid, p, fields)
    Tf = autotroph_temperature_scaling(a, T)

    return a.linear_mortality_rate * P′ * Tf
end

@inline function autotroph_aggregation(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)

    P′ = loss_active_biomass(val_prefix, i, j, k, grid, p, fields)

    return aggregation_rate(a.quadratic_mortality_rate, p.aggregation_exponent,
                            a.minimum_aggregation_rate, a.maximum_aggregation_rate, P′)
end

@inline autotroph_total_loss(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields) +
    autotroph_aggregation(val_prefix, i, j, k, grid, p, fields)

# the aggregate loss Σₐ = grazed + mortality + aggregation, which drives every biomass-element
# tendency through the variable quotas
@inline autotroph_aggregate_loss(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    zooplankton_grazing(val_prefix, i, j, k, grid, p, fields) +
    autotroph_total_loss(val_prefix, i, j, k, grid, p, fields)

#####
##### Uptake fluxes, calcification, oxygen production and loss routing
#####
##### These take the whole `bgc` because they need more than the plankton component.
#####

@inline photosynthesis_flux(val_prefix, i, j, k, grid, bgc, fields, aux) =
    photosynthesis(val_prefix, i, j, k, grid, bgc.plankton, fields, aux)

@inline photoacclimation_flux(val_prefix, i, j, k, grid, bgc, fields, aux) =
    photoacclimation(val_prefix, i, j, k, grid, bgc.plankton, fields, aux)

# N uptake, split between NO₃ and NH₄ by their Monod ratios. N-fixers take a total limitation of 1 (their
# N demand is met by fixation), which is a per-PFT trait ⇒ dispatched.
@inline function _nitrogen_monod(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    x = NO₃ / a.nutrient_half_saturations.nitrate
    y = max(zero(NH₄), NH₄ / a.nutrient_half_saturations.ammonia)  # floored, as in minimum_limitation

    return x, y, _total_nitrogen_limitation(Val(traits(a).nitrogen_fixer), x, y)
end

@inline _total_nitrogen_limitation(::Val{false}, x, y) = (x + y) / (1 + x + y)
@inline _total_nitrogen_limitation(::Val{true}, x, y)  = one(x)

@inline function nitrate_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    pc = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)

    x, y, VNtot = _nitrogen_monod(val_prefix, i, j, k, grid, p, fields)

    nitrate_fraction = x / (1 + x + y)

    return ifelse(VNtot > 0, nitrate_fraction / VNtot * pc * p.nitrogen_to_carbon, zero(pc))
end

@inline function ammonia_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    pc = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)

    x, y, VNtot = _nitrogen_monod(val_prefix, i, j, k, grid, p, fields)

    ammonia_fraction = y / (1 + x + y)

    return ifelse(VNtot > 0, ammonia_fraction / VNtot * pc * p.nitrogen_to_carbon, zero(pc))
end

# P uptake = photosynthesis·gQp, competitively partitioned between PO₄ and DOP by their Monod ratios. The
# autotroph P tracer sees the TOTAL; the split only routes which pool each part is drawn from.
@inline phosphorus_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux) =
    photosynthesis(val_prefix, i, j, k, grid, bgc.plankton, fields, aux) *
    growth_phosphorus_quota(val_prefix, i, j, k, grid, bgc.plankton, fields)

@inline function _phosphorus_monod(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    a = autotroph(p, val_prefix)

    @inbounds begin
        PO₄ = fields.PO₄[i, j, k]
        DOP = fields.DOP[i, j, k]
    end

    return PO₄ / a.nutrient_half_saturations.phosphate, DOP / a.nutrient_half_saturations.DOP
end

@inline function phosphate_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    u, v = _phosphorus_monod(val_prefix, i, j, k, grid, bgc.plankton, fields)

    Ptot = phosphorus_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)

    return ifelse(u + v > 0, u / (u + v) * Ptot, zero(Ptot))
end

@inline function dop_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    u, v = _phosphorus_monod(val_prefix, i, j, k, grid, bgc.plankton, fields)

    Ptot = phosphorus_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)

    return ifelse(u + v > 0, v / (u + v) * Ptot, zero(Ptot))
end

@inline iron_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux) =
    photosynthesis(val_prefix, i, j, k, grid, bgc.plankton, fields, aux) *
    growth_iron_quota(val_prefix, i, j, k, grid, bgc.plankton, fields)

@inline silicon_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux) =
    photosynthesis(val_prefix, i, j, k, grid, bgc.plankton, fields, aux) *
    growth_silicon_quota(val_prefix, i, j, k, grid, bgc.plankton, fields)

# the calcification FORMULA is a per-PFT trait: implicit (production fraction × photosynthesis × f_nut²,
# with a low-temperature taper and a bloom cap) or explicit (Krumhardt et al. 2019 picpoc)
@inline calcite_formation(val_prefix, i, j, k, grid, bgc, fields, aux) =
    _calcite_formation(Val(traits(bgc.plankton, val_prefix).explicit_calcifier),
                       val_prefix, i, j, k, grid, bgc, fields, aux)

@inline function _calcite_formation(::Val{false}, val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    @inbounds begin
        T = fields.T[i, j, k]
        C = autotroph_carbon(val_prefix, i, j, k, p, fields)
    end

    pc = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)
    fₙ = nutrient_limitation(val_prefix, i, j, k, grid, p, fields)

    return calcite_formation(p.calcite_production_fraction, p.maximum_calcifying_fraction,
                             p.calcite_bloom_threshold, p.calcite_temperature_threshold,
                             p.calcite_temperature_minimum, pc, fₙ, T, C)
end

@inline function _calcite_formation(::Val{true}, val_prefix, i, j, k, grid, bgc, fields, aux)
    p   = bgc.plankton
    a   = autotroph(p, val_prefix)
    khs = a.nutrient_half_saturations

    @inbounds begin
        T   = fields.T[i, j, k]
        PO₄ = fields.PO₄[i, j, k]
        DOP = fields.DOP[i, j, k]
        CO₂ = p.carbon_dioxide[i, j, k]
    end

    pc   = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)
    VPO₄ = phosphate_limitation(PO₄, DOP, khs.phosphate, khs.DOP)

    return picpoc_calcite_formation(CO₂, VPO₄, pc, T)
end

# N-fixation excretion to NH₄ (N-fixers only; zero otherwise, so it can be summed over all PFTs)
@inline nitrogen_excretion(val_prefix, i, j, k, grid, bgc, fields, aux) =
    _nitrogen_excretion(Val(traits(bgc.plankton, val_prefix).nitrogen_fixer),
                        val_prefix, i, j, k, grid, bgc, fields, aux)

@inline _nitrogen_excretion(::Val{false}, val_prefix, i, j, k, grid, bgc, fields, aux) =
    zero(bgc.plankton.nitrogen_to_carbon)

@inline function _nitrogen_excretion(::Val{true}, val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    pc    = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)
    NO₃_V = nitrate_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    NH₄_V = ammonia_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)

    return nitrogen_excretion(p.nitrogen_fixation_ratio, p.nitrogen_to_carbon, pc, NO₃_V, NH₄_V)
end

# O₂ produced by this PFT's carbon fixation, split by N source: nitrate-based new production releases
# more O₂ than ammonium-based regenerated production, and N-fixers add a fixation term
@inline function oxygen_production(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton
    o = bgc.oxygen

    pc    = photosynthesis(val_prefix, i, j, k, grid, p, fields, aux)
    NO₃_V = nitrate_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    NH₄_V = ammonia_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
    Nfix  = _nitrogen_fixation(Val(traits(p, val_prefix).nitrogen_fixer), p, pc, NO₃_V, NH₄_V)

    denominator = NO₃_V + NH₄_V + Nfix

    production = pc * (NO₃_V / denominator / o.nitrate_carbon_to_oxygen +
                       NH₄_V / denominator / o.remineralisation_carbon_to_oxygen +
                       Nfix  / denominator / o.diazotroph_carbon_to_oxygen)

    return ifelse(denominator > 0, production, zero(pc))
end

@inline _nitrogen_fixation(::Val{false}, p, pc, NO₃_V, NH₄_V) = zero(pc)

@inline _nitrogen_fixation(::Val{true}, p, pc, NO₃_V, NH₄_V) =
    nitrogen_fixation(p.nitrogen_fixation_ratio, p.nitrogen_to_carbon, pc, NO₃_V, NH₄_V)

#####
##### Loss routing
#####

# carbon: POC ← the particulate share of the linear loss plus aggregation; DIC/DOC split the remainder
@inline function autotroph_loss_to_dic(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    loss = autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields)

    return p.labile_fraction * (loss - loss_to_poc(val_prefix, i, j, k, grid, p, fields))
end

@inline function autotroph_loss_to_doc(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton

    loss = autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields)

    return (1 - p.labile_fraction) * (loss - loss_to_poc(val_prefix, i, j, k, grid, p, fields))
end

@inline autotroph_solid_carbon(val_prefix, i, j, k, grid, bgc, fields, aux) =
    loss_to_poc(val_prefix, i, j, k, grid, bgc.plankton, fields) +
    autotroph_aggregation(val_prefix, i, j, k, grid, bgc.plankton, fields)

# phosphorus: the P liberated by an autotroph's grazing/loss/aggregation (Qₚ·Σₐ) is split — the
# particulate share → POP, the grazed-to-zoo share is held by the zoo's fixed quota, and the surplus R
# splits to DOP and PO₄. If R < 0 the POP route is reduced.
@inline P_to_pop_gross(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    (grazing_to_poc(val_prefix, i, j, k, grid, p, fields) +
     loss_to_poc(val_prefix, i, j, k, grid, p, fields) +
     autotroph_aggregation(val_prefix, i, j, k, grid, p, fields)) *
    standing_phosphorus_quota(val_prefix, i, j, k, grid, p, fields)

@inline P_surplus(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    standing_phosphorus_quota(val_prefix, i, j, k, grid, p, fields) *
    autotroph_aggregate_loss(val_prefix, i, j, k, grid, p, fields) -
    p.zooplankton_phosphate_to_carbon *
    grazing_to_zooplankton(val_prefix, i, j, k, grid, p, fields) -
    P_to_pop_gross(val_prefix, i, j, k, grid, p, fields)

@inline function remaining_P_pop(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields)
    R = P_surplus(val_prefix, i, j, k, grid, p, fields)

    return P_to_pop_gross(val_prefix, i, j, k, grid, p, fields) + min(R, zero(R))
end

@inline remaining_P_dop(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    p.phosphorus_to_DOP_fraction *
    max(P_surplus(val_prefix, i, j, k, grid, p, fields), zero(p.nitrogen_to_carbon))

@inline remaining_P_dip(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    (1 - p.phosphorus_to_DOP_fraction) *
    max(P_surplus(val_prefix, i, j, k, grid, p, fields), zero(p.nitrogen_to_carbon))

# total biological Fe an autotroph liberates (i.e. not retained by the zooplankton)
@inline iron_loss(val_prefix, i, j, k, grid, bgc, fields, aux) =
    standing_iron_quota(val_prefix, i, j, k, grid, bgc.plankton, fields) *
    autotroph_aggregate_loss(val_prefix, i, j, k, grid, bgc.plankton, fields) -
    bgc.plankton.zooplankton_iron_to_carbon *
    grazing_to_zooplankton(val_prefix, i, j, k, grid, bgc.plankton, fields)

#####
##### The NutrientsPlanktonDetritus coupling contract, summed over autotrophs and zooplankton
#####

# dissolved-iron scavenging onto sinking particles belongs to the iron cycle, which lives in the
# `Nutrients` module — included before this one. The generic is declared here, where it is first needed,
# with a fallback of zero for a model with no iron cycle; `ComplexedIron` adds the real method.
@inline iron_scavenging(i, j, k, grid, iron, bgc, fields, aux) = zero(bgc.plankton.nitrogen_to_carbon)

# total primary production = gross carbon fixation summed over the PFTs. This drives the DIC drawdown
# directly, rather than overloading the carbon-currency `nutrient_uptake` base to mean photosynthesis.
@inline primary_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(photosynthesis_flux, i, j, k, grid, bgc, fields, aux)

@inline nutrient_uptake(i, j, k, grid, ::Val{:NO₃}, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(nitrate_uptake_flux, i, j, k, grid, bgc, fields, aux)

@inline nutrient_uptake(i, j, k, grid, ::Val{:NH₄}, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(ammonia_uptake_flux, i, j, k, grid, bgc, fields, aux)

@inline nutrient_uptake(i, j, k, grid, ::Val{:PO₄}, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(phosphate_uptake_flux, i, j, k, grid, bgc, fields, aux)

@inline nutrient_uptake(i, j, k, grid, ::Val{:Fe}, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(iron_uptake_flux, i, j, k, grid, bgc, fields, aux)

# carbon → DIC / DOC / POC: autotroph loss + autotroph grazing + zooplankton mortality AND predation
@inline inorganic_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(autotroph_loss_to_dic, i, j, k, grid, bgc, fields, aux) +
    sum_over_phytoplankton(grazing_to_dic, i, j, k, grid, p, fields) +
    sum_over_zooplankton(zooplankton_to_dic, i, j, k, grid, p, fields)

@inline dissolved_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(autotroph_loss_to_doc, i, j, k, grid, bgc, fields, aux) +
    sum_over_phytoplankton(grazing_to_doc, i, j, k, grid, p, fields) +
    sum_over_zooplankton(zooplankton_to_doc, i, j, k, grid, p, fields)

@inline solid_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(autotroph_solid_carbon, i, j, k, grid, bgc, fields, aux) +
    sum_over_phytoplankton(grazing_to_poc, i, j, k, grid, p, fields) +
    sum_over_zooplankton(zooplankton_to_poc, i, j, k, grid, p, fields)

# NH₄ source. N follows C at the fixed quota: of the biomass N liberated, the DIC-routed share
# mineralises straight to NH₄; of the DOC-routed share only the DON fraction is retained as DON (produced
# in the detritus), so the remainder mineralises to NH₄ here. Plus the N-fixers' excretion.
@inline inorganic_nitrogen_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    p.nitrogen_to_carbon * inorganic_waste(i, j, k, grid, p, bgc, fields, aux) +
    p.nitrogen_to_carbon * (1 - p.nitrogen_to_DON_fraction) *
        dissolved_waste(i, j, k, grid, p, bgc, fields, aux) +
    sum_over_nitrogen_fixers(nitrogen_excretion, i, j, k, grid, bgc, fields, aux)

# element-specific P/Fe waste to the INORGANIC pools. Variable quotas make the scalar-ratio decomposition
# invalid, so these cannot be derived from `phosphate_ratio`/`iron_ratio`. The autotrophs contribute only
# their dissolved-inorganic-routed surplus — the particulate and dissolved-organic shares go to the
# detritus POP/DOP tracers — and the zooplankton contribute at their fixed quota.
@inline inorganic_phosphate_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(remaining_P_dip, i, j, k, grid, p, fields) +
    p.zooplankton_phosphate_to_carbon *
        sum_over_zooplankton(zooplankton_to_dic, i, j, k, grid, p, fields)

# the biological Fe liberated (autotroph loss/aggregation/grazing not kept by a predator, plus the
# zooplankton's own liberated carbon at its quota) minus the particulate-routed share, which sinks as
# particulate iron instead
@inline inorganic_iron_waste(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(iron_loss, i, j, k, grid, bgc, fields, aux) +
    p.zooplankton_iron_to_carbon *
        sum_over_zooplankton(zooplankton_liberated_carbon, i, j, k, grid, p, fields) -
    pfe_production(i, j, k, grid, bgc, fields, aux)

# particulate-iron production. Because `inorganic_iron_waste` subtracts this and the PFe tendency adds it,
# folding scavenging in here makes ∂Fe pick up −scavenging and ∂PFe pick up +scavenging automatically, so
# Fe + PFe stays conserved.
@inline iron_to_pfe(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    standing_iron_quota(val_prefix, i, j, k, grid, p, fields) *
    (autotroph_aggregation(val_prefix, i, j, k, grid, p, fields) +
     grazing_to_poc(val_prefix, i, j, k, grid, p, fields) +
     loss_to_poc(val_prefix, i, j, k, grid, p, fields))

@inline pfe_production(i, j, k, grid, bgc, fields, aux) =
    sum_over_phytoplankton(iron_to_pfe, i, j, k, grid, bgc.plankton, fields) +
    bgc.plankton.zooplankton_iron_to_carbon *
        sum_over_zooplankton(zooplankton_to_poc, i, j, k, grid, bgc.plankton, fields) +
    iron_scavenging(i, j, k, grid, bgc.nutrients.iron, bgc, fields, aux)

#####
##### Scalar element:carbon ratios
#####
##### One value per component, but the quotas here are per-PFT and VARIABLE, so a scalar cannot carry the
##### plankton's P:C / Fe:C / Si:C — those flow through the element-specific waste functions above.
##### `carbon_ratio = 1` and `nitrogen_ratio = Q` ARE correct (carbon currency, fixed N:C).
#####

@inline carbon_ratio(i, j, k, grid, p::ManyPhytoZoo, ::NutrientsPlanktonDetritus, fields) =
    one(p.nitrogen_to_carbon)

@inline nitrogen_ratio(i, j, k, grid, p::ManyPhytoZoo, ::NutrientsPlanktonDetritus, fields) =
    p.nitrogen_to_carbon

@inline iron_ratio(i, j, k, grid, p::ManyPhytoZoo, ::NutrientsPlanktonDetritus, fields) =
    zero(p.nitrogen_to_carbon)

@inline silicon_ratio(i, j, k, grid, p::ManyPhytoZoo, ::NutrientsPlanktonDetritus, fields) =
    zero(p.nitrogen_to_carbon)

# decoupled from conservation, so its only live use is alkalinity's P:N. A default of zero omits the
# PO₄–alkalinity term; set `alkalinity_phosphate_ratio` (≈1/117) to include it, with no effect on C/P.
@inline phosphate_ratio(i, j, k, grid, p::ManyPhytoZoo, ::NutrientsPlanktonDetritus, fields) =
    p.alkalinity_phosphate_ratio

# total chlorophyll = the sum of the PFTs' Chl tracers, for the Chl-based light attenuation
@generated function chlorophyll(p::ManyPhytoZoo, model)
    exprs = [:(model.tracers.$(chlorophyll_name(s))) for s in phytoplankton_names(p)]
    return Expr(:call, :+, exprs...)
end

#####
##### Production terms consumed by a multi-element detritus
#####
##### These are sums over the plankton, but they feed the detritus tendencies. The generics are declared
##### in the parent module because the detritus is included first; the methods belong here, with the
##### machinery they reduce over.
#####

@inline dissolved_phosphorus_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    p.zooplankton_phosphate_to_carbon *
        sum_over_zooplankton(zooplankton_to_doc, i, j, k, grid, p, fields) +
    sum_over_phytoplankton(remaining_P_dop, i, j, k, grid, p, fields)

# the raw particulate P production may go negative; the shortfall is then taken from DOP instead
# (`dissolved_phosphorus_balance`) and this floored at zero, which keeps phosphorus conserved
@inline _particulate_phosphorus_production(i, j, k, grid, p::ManyPhytoZoo, fields) =
    p.zooplankton_phosphate_to_carbon *
        sum_over_zooplankton(zooplankton_to_poc, i, j, k, grid, p, fields) +
    sum_over_phytoplankton(remaining_P_pop, i, j, k, grid, p, fields)

@inline particulate_phosphorus_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    max(_particulate_phosphorus_production(i, j, k, grid, p, fields), zero(p.nitrogen_to_carbon))

@inline dissolved_phosphorus_balance(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    max(-_particulate_phosphorus_production(i, j, k, grid, p, fields), zero(p.nitrogen_to_carbon))

@inline dissolved_phosphorus_uptake(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(dop_uptake_flux, i, j, k, grid, bgc, fields, aux)

@inline particulate_iron_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    pfe_production(i, j, k, grid, bgc, fields, aux)

# silicate returned to solution: the grazed-and-dissolved share plus the non-particulate loss share, less
# this PFT's own uptake
@inline function _silicon_dissolution(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton
    a = autotroph(p, val_prefix)

    Qsi = standing_silicon_quota(val_prefix, i, j, k, grid, p, fields)

    return Qsi * (p.grazed_silica_remineralisation_fraction *
                      zooplankton_grazing(val_prefix, i, j, k, grid, p, fields) +
                  (1 - a.particulate_loss_fraction) *
                      autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields)) -
           silicon_uptake_flux(val_prefix, i, j, k, grid, bgc, fields, aux)
end

@inline silicon_dissolution(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_silicifiers(_silicon_dissolution, i, j, k, grid, bgc, fields, aux)

# the complementary share that sinks as biogenic silica
@inline function _biogenic_silica_production(val_prefix, i, j, k, grid, bgc, fields, aux)
    p = bgc.plankton
    a = autotroph(p, val_prefix)

    Qsi = standing_silicon_quota(val_prefix, i, j, k, grid, p, fields)

    return Qsi * ((1 - p.grazed_silica_remineralisation_fraction) *
                      zooplankton_grazing(val_prefix, i, j, k, grid, p, fields) +
                  autotroph_aggregation(val_prefix, i, j, k, grid, p, fields) +
                  a.particulate_loss_fraction *
                      autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields))
end

@inline biogenic_silica_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_silicifiers(_biogenic_silica_production, i, j, k, grid, bgc, fields, aux)

# total oxygen produced by carbon fixation, summed over the PFTs. Consumed by an oxygen component.
@inline oxygen_production_total(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_phytoplankton(oxygen_production, i, j, k, grid, bgc, fields, aux)

#####
##### Calcium carbonate
#####
##### The *living* calcite of each calcifying PFT is a plankton tracer (its tendency comes out of the
##### calcifier block of the manifest); the *detrital* calcite is the `ExplicitCalciumCarbonate` pool in
##### the inorganic-carbon slot. Rather than replacing that component's whole net production, we override
##### its three per-plankton hooks, each summed over the calcifying PFTs, so it keeps assembling DIC and
##### alkalinity itself. Conservation then closes across the living pool too:
##### ∂(living CaCO₃) + ∂CaCO₃ + ∂DIC_calcite = 0, with loss, aggregation and grazing all cancelling.
#####

const ManyPhytoZoo_NPD_EC =
    NutrientsPlanktonDetritus{<:Any, <:Nutrients, <:ManyPhytoZoo, <:Any, <:ExplicitCalciumCarbonate}

# calcite routed to the sinking pool: the quota times the losses plus the share of grazed calcite that is
# not dissolved in the gut
@inline _calcite_to_pool(val_prefix, i, j, k, grid, bgc, fields, aux) =
    standing_calcite_quota(val_prefix, i, j, k, grid, bgc.plankton, fields) *
    (autotroph_total_loss(val_prefix, i, j, k, grid, bgc.plankton, fields) +
     (1 - bgc.plankton.grazed_calcite_remineralisation_fraction) *
         zooplankton_grazing(val_prefix, i, j, k, grid, bgc.plankton, fields))

# grazed calcite dissolved in the gut and returned to DIC and alkalinity
@inline _calcite_gut_dissolution(val_prefix, i, j, k, grid, bgc, fields, aux) =
    bgc.plankton.grazed_calcite_remineralisation_fraction *
    standing_calcite_quota(val_prefix, i, j, k, grid, bgc.plankton, fields) *
    zooplankton_grazing(val_prefix, i, j, k, grid, bgc.plankton, fields)

@inline biological_calcium_carbonate_precipitation(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_calcifiers(calcite_formation, i, j, k, grid, bgc, fields, aux)

@inline particulate_calcium_carbonate_production(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_calcifiers(_calcite_to_pool, i, j, k, grid, bgc, fields, aux)

@inline biological_calcium_carbonate_dissolution(i, j, k, grid, p::ManyPhytoZoo, bgc, fields, aux) =
    sum_over_calcifiers(_calcite_gut_dissolution, i, j, k, grid, bgc, fields, aux)

# a first-order pool dissolution with no saturation term, overriding the component's saturation-law
# default for this plankton. The rate lives on `ExplicitCalciumCarbonate`.
@inline function calcium_carbonate_dissolution_flux(i, j, k, grid, bgc::ManyPhytoZoo_NPD_EC,
                                                    fields, auxiliary_fields,
                                                    ::Val{calcite_name}, ::Val{Ω_name},
                                                    ::Val{n}) where {calcite_name, Ω_name, n}
    @inbounds CaCO₃ = getproperty(fields, calcite_name)[i, j, k]

    return @inbounds bgc.inorganic_carbon.calcium_carbonate_dissolution_rate[n] * CaCO₃
end

#####
##### Aqueous CO₂
#####
##### A field owned by the plankton and recomputed each step from the PRIMARY dissolved-inorganic-carbon
##### and alkalinity realization through the inorganic-carbon solver. The carbon-limited and explicit-
##### calcifier rate laws read it directly off the plankton. Reading only the primary realization keeps
##### the plankton↔carbonate coupling single-realization when several are carried.
#####

update_biogeochemical_state!(model, p::ManyPhytoZoo, npd::NutrientsPlanktonDetritus) =
    update_carbon_dioxide!(model, p.carbon_dioxide, npd)

# no carbon-limited or explicit-calcifier PFT ⇒ no field to fill
update_carbon_dioxide!(model, ::Nothing, npd) = nothing

function update_carbon_dioxide!(model, carbon_dioxide, npd)
    grid = model.grid

    DIC_name, Alk_name = carbonate_field_names(npd.inorganic_carbon)[1]

    launch!(architecture(grid), grid, :xyz, _compute_carbon_dioxide!,
            carbon_dioxide, npd.inorganic_carbon.carbon_chemistry, grid, fields(model),
            Val(DIC_name), Val(Alk_name), defaults.gravitational_acceleration)

    fill_halo_regions!(carbon_dioxide, grid, model.clock, fields(model))

    return nothing
end

@kernel function _compute_carbon_dioxide!(carbon_dioxide, carbon_chemistry, grid, model_fields,
                                          ::Val{DIC_name}, ::Val{Alk_name}, g) where {DIC_name, Alk_name}
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        T   = model_fields.T[i, j, k]
        S   = model_fields.S[i, j, k]
        DIC = getproperty(model_fields, DIC_name)[i, j, k]
        Alk = getproperty(model_fields, Alk_name)[i, j, k]

        # total alkalinity enters the pH solve, so phosphate and silicate are needed too
        phosphate = model_fields.PO₄[i, j, k]
        silicate  = model_fields.Si[i, j, k]
    end

    P = hydrostatic_pressure(model_fields, i, j, k, grid, g)

    CO₂ = carbon_chemistry(; DIC, Alk, T, S, P, phosphate, silicate, output = Val(:CO₂))

    @inbounds carbon_dioxide[i, j, k] = CO₂
end
