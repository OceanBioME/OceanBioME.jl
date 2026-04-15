@kwdef struct NanoAndDiatoms{N, D, FT}
               nano :: N
            diatoms :: D
    base_rain_ratio :: FT = 0.3
end

required_biogeochemical_tracers(phyto::NanoAndDiatoms) = (required_biogeochemical_tracers(phyto.nano, :P)...,
                                                          required_biogeochemical_tracers(phyto.diatoms, :D)...)

@inline dissolved_exudate(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (dissolved_exudate(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + dissolved_exudate(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline function uptake(nano_exudated_fraction,
                       nano_maximum_iron_ratio,
                       nano_half_saturation_for_iron_uptake,
                       nano_threshold_for_size_dependency,
                       nano_size_ratio,
                       nano_minimum_ammonium_half_saturation,
                       nano_minimum_nitrate_half_saturation,
                       nano_optimal_iron_quota,
                       nano_base_growth_rate,
                       nano_temperature_sensitivity,
                       diatom_exudated_fraction,
                       diatom_maximum_iron_ratio,
                       diatom_half_saturation_for_iron_uptake,
                       diatom_threshold_for_size_dependency,
                       diatom_size_ratio,
                       diatom_minimum_ammonium_half_saturation,
                       diatom_minimum_nitrate_half_saturation,
                       diatom_optimal_iron_quota,
                       diatom_base_growth_rate,
                       diatom_temperature_sensitivity,
                       ::Val{:Fe},
                       T,
                       Fe,
                       NO₃,
                       NH₄,
                       PO₄,
                       Si,
                       Si′,
                       P,
                       PChl,
                       PFe,
                       D,
                       DChl,
                       DFe)
    return (iron_uptake(nano_exudated_fraction,
                        nano_maximum_iron_ratio,
                        nano_half_saturation_for_iron_uptake,
                        nano_threshold_for_size_dependency,
                        nano_size_ratio,
                        nano_minimum_ammonium_half_saturation,
                        nano_minimum_nitrate_half_saturation,
                        nano_optimal_iron_quota,
                        nano_base_growth_rate,
                        nano_temperature_sensitivity,
                        T,
                        Fe,
                        P,
                        PChl,
                        PFe,
                        NO₃,
                        NH₄)
            + iron_uptake(diatom_exudated_fraction,
                          diatom_maximum_iron_ratio,
                          diatom_half_saturation_for_iron_uptake,
                          diatom_threshold_for_size_dependency,
                          diatom_size_ratio,
                          diatom_minimum_ammonium_half_saturation,
                          diatom_minimum_nitrate_half_saturation,
                          diatom_optimal_iron_quota,
                          diatom_base_growth_rate,
                          diatom_temperature_sensitivity,
                          T,
                          Fe,
                          D,
                          DChl,
                          DFe,
                          NO₃,
                          NH₄))
end

@inline uptake(phyto::NanoAndDiatoms,
               ::Val{:Fe},
               T,
               Fe,
               NO₃,
               NH₄,
               PO₄,
               Si,
               Si′,
               P,
               PChl,
               PFe,
               D,
               DChl,
               DFe) =
    uptake(phyto.nano.exudated_fraction,
           phyto.nano.maximum_iron_ratio,
           phyto.nano.half_saturation_for_iron_uptake,
           phyto.nano.threshold_for_size_dependency,
           phyto.nano.size_ratio,
           phyto.nano.nutrient_limitation.minimum_ammonium_half_saturation,
           phyto.nano.nutrient_limitation.minimum_nitrate_half_saturation,
           phyto.nano.nutrient_limitation.optimal_iron_quota,
           phyto.nano.growth_rate.base_growth_rate,
           phyto.nano.growth_rate.temperature_sensitivity,
           phyto.diatoms.exudated_fraction,
           phyto.diatoms.maximum_iron_ratio,
           phyto.diatoms.half_saturation_for_iron_uptake,
           phyto.diatoms.threshold_for_size_dependency,
           phyto.diatoms.size_ratio,
           phyto.diatoms.nutrient_limitation.minimum_ammonium_half_saturation,
           phyto.diatoms.nutrient_limitation.minimum_nitrate_half_saturation,
           phyto.diatoms.nutrient_limitation.optimal_iron_quota,
           phyto.diatoms.growth_rate.base_growth_rate,
           phyto.diatoms.growth_rate.temperature_sensitivity,
           Val(:Fe),
           T,
           Fe,
           NO₃,
           NH₄,
           PO₄,
           Si,
           Si′,
           P,
           PChl,
           PFe,
           D,
           DChl,
           DFe)

@inline uptake(phyto::NanoAndDiatoms, val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (uptake(phyto.nano, Val(:P), val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + uptake(phyto.diatoms, Val(:D), val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline function nitrogen_availability_limitation(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    _, _, _, LN = phyto.nano.nutrient_limitation(Val(:P), i, j, k, grid, bgc, phyto.nano, clock, fields, auxiliary_fields)

    return LN
end

@inline base_production_rate(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    @inbounds base_production_rate(phyto.nano.growth_rate, fields.T[i, j, k])

@inline silicate_uptake(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    silicate_uptake(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline total_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (total_production(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + total_production(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields))