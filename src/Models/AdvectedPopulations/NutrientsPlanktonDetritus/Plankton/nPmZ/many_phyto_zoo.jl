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

const _manifested_marbl_plankton = Set{Tuple}()

@inline autotroph(p::Autotroph, ::Val{S}) where S = getproperty(p.autotrophs, S)
@inline traits(p::Autotroph, val_sname) = traits(autotroph(p, val_sname))


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

@inline function nutrient_limitation(val_prefix, i, j, k, grid, p::MARBLPlankton, fields)
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

    f = minimum_limitation(traits(a).silicifier, traits(a).nfixer,
                           khs.nitrate, khs.ammonia, khs.phosphate, khs.DOP, khs.iron, khs.silicate,
                           NO₃, NH₄, PO₄, DOP, Fe, SiO₃)

    return carbon_limitation(Val(traits(a).carbon_limited), f, a, p, i, j, k)
end