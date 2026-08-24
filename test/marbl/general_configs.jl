# Generated from marbl-info3/draft/general_config/configs.jl — the general N-plankton configurations
# shared by the comparison harnesses and the MARBL settings files that produced their baselines.

#####
##### Phase 10 general-config definitions — the SINGLE SOURCE OF TRUTH for both sides of the comparison.
#####
##### These two configs exist to exercise the general N_A/N_Z machinery, and each is used twice: to build
##### the OceanBioME model (`general2_plankton` / `general5_plankton`) AND to generate the MARBL Fortran
##### settings file (write_marbl_settings.jl reads the SAME structs), so the two implementations cannot
##### drift apart by a typo.
#####
##### `general5` is the demanding one — it turns on every code path Phase 10 adds:
#####   • 5 autotrophs: TWO calcifiers (sp implicit §9.1, cocco explicit §9.2) and TWO silicifiers
#####     (diat, diat2) ⇒ sum_over_calcifiers / sum_over_silicifiers really reduce over >1 PFT
#####   • 2 zooplankton (zoo, zoo2) ⇒ x_graze_zoo, f_zoo_detr and the zoo tendencies are per-predator
#####   • MULTI-MEMBER prey classes: zoo eats (sp, cocco) and (diaz, diat2) as single pools whose grazing
#####     flux is split by biomass share (MARBL auto_ind_cnt = 2)
#####   • ZOO-ON-ZOO grazing: zoo2 eats zoo (MARBL zoo_ind_cnt = 1) ⇒ the zoo_graze_* terms, which are
#####     identically zero in every CESM config and so were never implemented before Phase 10
#####   • prey grazed by BOTH predators (sp, diat2) ⇒ auto_graze/auto_graze_zootot sum over the matrix
#####
##### `general2` is the degenerate end: 2 autotrophs, 1 zooplankton, NO N-fixer and NO carbon-limited PFT
##### (so sum_over_nfixers reduces over the empty set and `carbon_dioxide` is `nothing`).
#####
##### Parameters are CESM2.1 for sp/diat/diaz/zoo and the CESM2.1+cocco set for cocco; diat2/zoo2 are
##### deliberately-perturbed variants (a config that gave two PFTs identical parameters could not tell a
##### per-PFT bug from a shared one). The numbers themselves are arbitrary — what the gate tests is that
##### OceanBioME and MARBL agree given the SAME numbers.
#####

using OceanBioME, Oceananigans
using Oceananigans.Fields: CenterField

# --- the 5-autotroph / 2-zooplankton config -------------------------------------------------------

# a second silicifier: a slower, more iron-limited diatom
general_diatoms_2(FT = Float64) = MARBLAutotroph(FT;
    silicifier = true,
    nitrate = 0.8, ammonia = 0.08, phosphate = 0.06, DOP = 0.6, iron = 9.0e-5, silicate = 1.1,
    photosynthesis_rate_reference_per_day = 3.6, initial_PI_slope_per_day = 0.32,
    maximum_chlorophyll_to_nitrogen = 3.8, temperature_threshold = -10.0,
    quadratic_mortality_rate_per_day = 0.012, minimum_aggregation_rate_per_day = 0.02,
    loss_threshold_concentration = 0.025, cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 40.0e-6, minimum_iron_quota = 3.0e-6)

# predator 1: eats (sp, cocco) as one prey class, diat alone, and (diaz, diat2) as one class
function general_zooplankton_1(FT = Float64)
    grazing = (
        sp_cocco   = MARBLGrazing(FT; prey = (:sp, :cocco), maximum_grazing_rate_per_day = 3.3,
                                      grazing_half_saturation = 1.2, fraction_to_zooplankton = 0.3,
                                      fraction_to_particulate = 0.0, detritus_fraction = 0.12),
        diat       = MARBLGrazing(FT; prey = :diat, maximum_grazing_rate_per_day = 3.15,
                                      grazing_half_saturation = 1.2, fraction_to_zooplankton = 0.25,
                                      fraction_to_particulate = 0.39, detritus_fraction = 0.24),
        diaz_diat2 = MARBLGrazing(FT; prey = (:diaz, :diat2), maximum_grazing_rate_per_day = 3.3,
                                      grazing_half_saturation = 1.05, fraction_to_zooplankton = 0.3,
                                      fraction_to_particulate = 0.1, detritus_fraction = 0.12))

    return MARBLZooplankton(FT; grazing)
end

# predator 2: eats diat2, the OTHER ZOOPLANKTON (zoo), and sp — sigmoidal on its zooplankton prey
function general_zooplankton_2(FT = Float64)
    grazing = (
        diat2 = MARBLGrazing(FT; prey = :diat2, maximum_grazing_rate_per_day = 2.4,
                                 grazing_half_saturation = 0.9, fraction_to_zooplankton = 0.28,
                                 fraction_to_particulate = 0.35, detritus_fraction = 0.2),
        zoo   = MARBLGrazing(FT; prey = :zoo, maximum_grazing_rate_per_day = 1.8,
                                 grazing_half_saturation = 0.75, fraction_to_zooplankton = 0.32,
                                 fraction_to_particulate = 0.25, detritus_fraction = 0.15,
                                 sigmoidal = true),
        sp    = MARBLGrazing(FT; prey = :sp, maximum_grazing_rate_per_day = 2.1,
                                 grazing_half_saturation = 1.1, fraction_to_zooplankton = 0.26,
                                 fraction_to_particulate = 0.05, detritus_fraction = 0.1))

    # NB the Q10 (1.7) and its reference temperature (30 °C) are NOT perturbed here, and cannot be: MARBL
    # holds Q_10 as a module constant (marbl_settings_mod L97) and DERIVES Tref from temp_func_form_opt
    # (set_derived_from_temp_func_form), so neither is a per-PFT setting. Only the mortality/threshold
    # parameters differ from `zoo`.
    return MARBLZooplankton(FT; grazing,
                            linear_mortality_rate_per_day    = 0.12,
                            quadratic_mortality_rate_per_day = 0.35,
                            loss_threshold_concentration     = 0.06)
end

# the autotroph ORDER here is MARBL's auto_ind order (1..5); the zooplankton order is its zoo_ind order
general5_autotrophs(FT = Float64) = (sp    = default_small_phytoplankton(FT),
                                     diat  = default_diatoms(FT),
                                     diaz  = default_diazotrophs(FT),
                                     cocco = default_coccolithophores(FT),
                                     diat2 = general_diatoms_2(FT))

general5_zooplankton(FT = Float64) = (zoo  = general_zooplankton_1(FT),
                                      zoo2 = general_zooplankton_2(FT))

# needs a grid: cocco is carbon-limited, so the plankton holds the aqueous-CO₂ field its VCO2/picpoc read
general5_plankton(grid, FT = Float64) =
    MARBLPlankton(FT; autotrophs = general5_autotrophs(FT), zooplankton = general5_zooplankton(FT),
                      carbon_dioxide = CenterField(grid))

# --- the 2-autotroph / 1-zooplankton config -------------------------------------------------------

general2_autotrophs(FT = Float64) = (sp   = default_small_phytoplankton(FT),
                                     diat = default_diatoms(FT))

function general2_zooplankton(FT = Float64)
    grazing = (
        sp   = MARBLGrazing(FT; prey = :sp,   maximum_grazing_rate_per_day = 3.3,
                                fraction_to_zooplankton = 0.3, fraction_to_particulate = 0.0,
                                detritus_fraction = 0.12),
        diat = MARBLGrazing(FT; prey = :diat, maximum_grazing_rate_per_day = 3.15,
                                fraction_to_zooplankton = 0.25, fraction_to_particulate = 0.39,
                                detritus_fraction = 0.24))

    return (zoo = MARBLZooplankton(FT; grazing), )
end

general2_plankton(FT = Float64) =
    MARBLPlankton(FT; autotrophs = general2_autotrophs(FT), zooplankton = general2_zooplankton(FT))
