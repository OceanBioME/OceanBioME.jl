"""
Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model.

Tracers
=======

* Nano-phytoplankton: P (μmolC/L)
* Diatoms: D (μmolC/L)
* Zooplankton: Z (μmolC/L)
* Mesozooplankton: M (μmolC/L)
* Chlorophyll in nano-phytoplankton: Pᶜʰˡ (μgChl/L)
* Chlorophyll in diatoms: Dᶜʰˡ (μgChl/L)
* Iron in nano-phytoplanktons: Pᶠᵉ (nmolFe/L) 
* Iron in diatoms: Dᶠᵉ (nmolFe/L) 
* Silicon in diatoms: Dˢⁱ (μmolSi/L)

* Dissolved organic carbon: DOC (μmolC/L)
* Small sinking particles : POC (μmolC/L)
* Large sinking particles: GOC (μmolC/L)
* Iron in small particles: SFe (nmolFe/L) 
* Iron in large particles: BFe (nmolFe/L) 
* Silicate in large particles : PSi (μmolSi/L)
* Nitrates: NO₃ (μmolN/L)
* Ammonium: NH₄ (μmolN/L)
* Phosphate: PO₄ (μmolP/L)
* Dissolved iron: Fe (nmolFe/L)
* Silicate: Si (μmolSi/L)
* Calcite: CaCO₃ (μmolC/L)
* Dissolved oxygen: O₂ (μmolO₂/L)
* Dissolved inorganic carbon: DIC (μmolC/L)
* Total alkalinity: Alk (μmolN/L)


Required submodels
==================
# you will need something like this, they use a different PAR model but I wouldn't worry about that for now, you might also need temperatre and salinity (not sure)
* Photosynthetically available radiation: PAR (W/m²)

"""
module PISCESModel

export PISCES

using Oceananigans.Units
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField, ConstantField

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR
using OceanBioME: setup_velocity_fields, show_sinking_velocities, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.BoxModels: BoxModel
using OceanBioME.Boundaries.Sediments: sinking_flux

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

import OceanBioME: redfield, conserved_tracers

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity,
                                     biogeochemical_auxiliary_fields

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt
import Base: show, summary

import OceanBioME.Boundaries.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

struct PISCES{FT, PD, ZM, OT, W, CF, ZF} <: AbstractContinuousFormBiogeochemistry

    growth_rate_at_zero :: FT # add list of parameters here, assuming theyre all just numbers FT will be fine for advect_particles_kernel
    growth_rate_reference_for_light_limitation :: FT
    basal_respiration_rate :: FT
    temperature_sensitivity_of_growth :: FT
    initial_slope_of_PI_curve :: PD
    exudation_of_DOC :: PD
    absorption_in_the_blue_part_of_light :: PD
    absorption_in_the_green_part_of_light :: PD
    absorption_in_the_red_part_of_light :: PD
    min_half_saturation_const_for_phosphate :: PD
    min_half_saturation_const_for_ammonium :: PD
    min_half_saturation_const_for_nitrate :: PD
    min_half_saturation_const_for_silicate :: FT #
    parameter_for_half_saturation_const :: FT
    parameter_for_SiC :: OT #
    min_half_saturation_const_for_iron_uptake :: PD #
    size_ratio_of_phytoplankton :: PD#
    optimal_SiC_uptake_ratio_of_diatoms :: FT
    optimal_iron_quota :: PD#
    max_iron_quota :: PD#
    phytoplankton_mortality_rate :: PD
    min_quadratic_mortality_of_phytoplankton :: FT
    max_quadratic_mortality_of_diatoms :: FT
    max_ChlC_ratios_of_phytoplankton :: PD
    min_ChlC_ratios_of_phytoplankton :: FT
    threshold_concentration_for_size_dependency :: PD
    mean_residence_time_of_phytoplankton_in_unlit_mixed_layer :: PD

    latitude :: FT
    length_of_day :: FT

    temperature_sensitivity_term :: ZM
    max_growth_efficiency_of_zooplankton :: ZM
    non_assimilated_fraction :: ZM
    excretion_as_DOM :: ZM
    max_grazing_rate :: ZM
    flux_feeding_rate :: FT
    half_saturation_const_for_grazing :: ZM
    preference_for_nanophytoplankton :: ZM
    preference_for_diatoms :: ZM
    preference_for_POC :: ZM
    preference_for_microzooplankton :: FT
    food_threshold_for_zooplankton :: ZM
    specific_food_thresholds_for_microzooplankton :: FT
    specific_food_thresholds_for_mesozooplankton :: FT
    zooplankton_quadratic_mortality :: ZM
    zooplankton_linear_mortality :: ZM
    half_saturation_const_for_mortality :: FT
    fraction_of_calcite_not_dissolving_in_guts :: ZM
    FeC_ratio_of_zooplankton :: FT
    FeZ_redfield_ratio :: FT


    remineralisation_rate_of_DOC :: FT
    half_saturation_const_for_DOC_remin :: FT
    NO3_half_saturation_const_for_DOC_remin :: FT
    NH4_half_saturation_const_for_DOC_remin :: FT
    PO4_half_saturation_const_for_DOC_remin :: FT
    Fe_half_saturation_const_for_DOC_remin :: FT
    aggregation_rate_of_DOC_to_POC_1 :: FT
    aggregation_rate_of_DOC_to_POC_2 :: FT
    aggregation_rate_of_DOC_to_GOC_3 :: FT
    aggregation_rate_of_DOC_to_POC_4 :: FT
    aggregation_rate_of_DOC_to_POC_5 :: FT


    degradation_rate_of_POC :: FT
    sinking_speed_of_POC :: FT
    min_sinking_speed_of_GOC :: FT
    sinking_speed_of_dust :: FT
    aggregation_rate_of_POC_to_GOC_6 :: FT
    aggregation_rate_of_POC_to_GOC_7 :: FT
    aggregation_rate_of_POC_to_GOC_8 :: FT
    aggregation_rate_of_POC_to_GOC_9 :: FT
    min_scavenging_rate_of_iron :: FT
    slope_of_scavenging_rate_of_iron :: FT
    scavenging_rate_of_iron_by_dust :: FT
    dissolution_rate_of_calcite :: FT
    exponent_in_the_dissolution_rate_of_calcite :: FT
    proportion_of_the_most_labile_phase_in_PSi :: FT
    slow_dissolution_rate_of_PSi :: FT
    fast_dissolution_rate_of_PSi :: FT


    max_nitrification_rate :: FT
    half_sat_const_for_denitrification1 :: FT
    half_sat_const_for_denitrification2 :: FT
    total_concentration_of_iron_ligands :: FT
    max_rate_of_nitrogen_fixation :: FT
    Fe_half_saturation_constant_of_nitrogen_fixation :: FT
    photosynthetic_parameter_of_nitrogen_fixation :: FT
    iron_concentration_in_sea_ice :: FT
    max_sediment_flux_of_Fe :: FT
    solubility_of_iron_in_dust :: FT
    OC_for_ammonium_based_processes :: FT
    OC_ratio_of_nitrification :: FT
    CN_ratio_of_ammonification :: FT
    CN_ratio_of_denitrification :: FT
    NC_redfield_ratio :: FT
    PC_redfield_ratio :: FT
    rain_ratio_parameter :: FT
    bacterial_reference :: FT

    NC_stoichiometric_ratio_of_dentitrification :: FT
    NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER :: FT
    dissolution_rate_of_silicon :: FT
    coefficient_of_bacterial_uptake_of_iron_in_POC :: FT
    coefficient_of_bacterial_uptake_of_iron_in_GOC :: FT
    max_FeC_ratio_of_bacteria :: FT
    Fe_half_saturation_const_for_PLACEHOLDER :: FT    #not sure what this should be called
    proportion_of_sinking_grazed_shells :: ZM

    mixed_layer_depth :: CF
    euphotic_layer_depth :: CF
    yearly_maximum_silicate :: FT
    dust_deposition :: FT

    vertical_diffusivity :: CF 
    carbonate_sat_ratio :: ZF

    sinking_velocities :: W

    function PISCES(growth_rate_at_zero :: FT,
                    growth_rate_reference_for_light_limitation :: FT,
                    basal_respiration_rate :: FT,
                    temperature_sensitivity_of_growth :: FT,
                    initial_slope_of_PI_curve :: PD,
                    exudation_of_DOC :: PD,
                    absorption_in_the_blue_part_of_light :: PD,
                    absorption_in_the_green_part_of_light :: PD,
                    absorption_in_the_red_part_of_light :: PD,
                    min_half_saturation_const_for_phosphate :: PD,
                    min_half_saturation_const_for_ammonium :: PD,
                    min_half_saturation_const_for_nitrate :: PD,
                    min_half_saturation_const_for_silicate :: FT,
                    parameter_for_half_saturation_const :: FT,
                    parameter_for_SiC :: OT,
                    min_half_saturation_const_for_iron_uptake :: PD,
                    size_ratio_of_phytoplankton :: PD,
                    optimal_SiC_uptake_ratio_of_diatoms :: FT,
                    optimal_iron_quota :: PD,
                    max_iron_quota :: PD,
                    phytoplankton_mortality_rate :: PD,
                    min_quadratic_mortality_of_phytoplankton :: FT,
                    max_quadratic_mortality_of_diatoms :: FT,
                    max_ChlC_ratios_of_phytoplankton :: PD,
                    min_ChlC_ratios_of_phytoplankton :: FT,
                    threshold_concentration_for_size_dependency :: PD,
                    mean_residence_time_of_phytoplankton_in_unlit_mixed_layer :: PD,

                    latitude :: FT,
                    length_of_day :: FT, 
        
                    temperature_sensitivity_term :: ZM,
                    max_growth_efficiency_of_zooplankton :: ZM,
                    non_assimilated_fraction :: ZM,
                    excretion_as_DOM :: ZM,
                    max_grazing_rate :: ZM,
                    flux_feeding_rate :: FT,
                    half_saturation_const_for_grazing :: ZM,
                    preference_for_nanophytoplankton :: ZM,
                    preference_for_diatoms :: ZM,
                    preference_for_POC :: ZM,
                    preference_for_microzooplankton :: FT,
                    food_threshold_for_zooplankton :: ZM,
                    specific_food_thresholds_for_microzooplankton :: FT,
                    specific_food_thresholds_for_mesozooplankton :: FT,
                    zooplankton_quadratic_mortality :: ZM,
                    zooplankton_linear_mortality :: ZM,
                    half_saturation_const_for_mortality :: FT,
                    fraction_of_calcite_not_dissolving_in_guts :: ZM,
                    FeC_ratio_of_zooplankton :: FT,
                    FeZ_redfield_ratio :: FT,
    
    
                    remineralisation_rate_of_DOC :: FT,
                    half_saturation_const_for_DOC_remin :: FT,
                    NO3_half_saturation_const_for_DOC_remin :: FT,
                    NH4_half_saturation_const_for_DOC_remin :: FT,
                    PO4_half_saturation_const_for_DOC_remin :: FT,
                    Fe_half_saturation_const_for_DOC_remin :: FT,
                    aggregation_rate_of_DOC_to_POC_1 :: FT,
                    aggregation_rate_of_DOC_to_POC_2 :: FT,
                    aggregation_rate_of_DOC_to_GOC_3 :: FT,
                    aggregation_rate_of_DOC_to_POC_4 :: FT,
                    aggregation_rate_of_DOC_to_POC_5 :: FT,
    
    
                    degradation_rate_of_POC :: FT,
                    sinking_speed_of_POC :: FT,
                    min_sinking_speed_of_GOC :: FT,
                    sinking_speed_of_dust :: FT,
                    aggregation_rate_of_POC_to_GOC_6 :: FT,
                    aggregation_rate_of_POC_to_GOC_7 :: FT,
                    aggregation_rate_of_POC_to_GOC_8 :: FT,
                    aggregation_rate_of_POC_to_GOC_9 :: FT,
                    min_scavenging_rate_of_iron :: FT,
                    slope_of_scavenging_rate_of_iron :: FT,
                    scavenging_rate_of_iron_by_dust :: FT,
                    dissolution_rate_of_calcite :: FT,
                    exponent_in_the_dissolution_rate_of_calcite :: FT,
                    proportion_of_the_most_labile_phase_in_PSi :: FT,
                    slow_dissolution_rate_of_PSi :: FT,
                    fast_dissolution_rate_of_PSi :: FT,
    
    
                    max_nitrification_rate :: FT,
                    half_sat_const_for_denitrification1 :: FT,
                    half_sat_const_for_denitrification2 :: FT,
                    total_concentration_of_iron_ligands :: FT,
                    max_rate_of_nitrogen_fixation :: FT,
                    Fe_half_saturation_constant_of_nitrogen_fixation :: FT,
                    photosynthetic_parameter_of_nitrogen_fixation :: FT,
                    iron_concentration_in_sea_ice :: FT,
                    max_sediment_flux_of_Fe :: FT,
                    solubility_of_iron_in_dust :: FT,
                    OC_for_ammonium_based_processes :: FT,
                    OC_ratio_of_nitrification :: FT,
                    CN_ratio_of_ammonification :: FT,
                    CN_ratio_of_denitrification :: FT,
                    NC_redfield_ratio :: FT,
                    PC_redfield_ratio :: FT,
                    rain_ratio_parameter :: FT,
                    bacterial_reference :: FT, 

                    NC_stoichiometric_ratio_of_dentitrification :: FT,
                    NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER :: FT,
                    dissolution_rate_of_silicon :: FT,
                    coefficient_of_bacterial_uptake_of_iron_in_POC :: FT,
                    coefficient_of_bacterial_uptake_of_iron_in_GOC :: FT,
                    max_FeC_ratio_of_bacteria :: FT,
                    Fe_half_saturation_const_for_PLACEHOLDER :: FT,    #not sure what this should be called
                    proportion_of_sinking_grazed_shells :: ZM,
                    
                    mixed_layer_depth :: CF,
                    euphotic_layer_depth :: CF,
                    yearly_maximum_silicate :: FT,
                    dust_deposition :: FT,
                    vertical_diffusivity :: CF, 
                    carbonate_sat_ratio :: ZF,

                    sinking_velocities :: W,) where {FT, PD, ZM, OT, W, CF, ZF} # then do the same here (this is all just annoying boiler plate but we need it to make the next function work)


        return new{FT, PD, ZM, OT, W, CF, ZF}(growth_rate_at_zero,
                            growth_rate_reference_for_light_limitation,
                            basal_respiration_rate,
                            temperature_sensitivity_of_growth,
                            initial_slope_of_PI_curve,
                            exudation_of_DOC,
                            absorption_in_the_blue_part_of_light,
                            absorption_in_the_green_part_of_light,
                            absorption_in_the_red_part_of_light,
                            min_half_saturation_const_for_phosphate,
                            min_half_saturation_const_for_ammonium,
                            min_half_saturation_const_for_nitrate,
                            min_half_saturation_const_for_silicate,
                            parameter_for_half_saturation_const,
                            parameter_for_SiC,
                            min_half_saturation_const_for_iron_uptake,
                            size_ratio_of_phytoplankton,
                            optimal_SiC_uptake_ratio_of_diatoms,
                            optimal_iron_quota,
                            max_iron_quota,
                            phytoplankton_mortality_rate,
                            min_quadratic_mortality_of_phytoplankton,
                            max_quadratic_mortality_of_diatoms,
                            max_ChlC_ratios_of_phytoplankton,
                            min_ChlC_ratios_of_phytoplankton,
                            threshold_concentration_for_size_dependency,
                            mean_residence_time_of_phytoplankton_in_unlit_mixed_layer,

                            latitude, 
                            length_of_day,

                            temperature_sensitivity_term,
                            max_growth_efficiency_of_zooplankton,
                            non_assimilated_fraction,
                            excretion_as_DOM,
                            max_grazing_rate,
                            flux_feeding_rate,
                            half_saturation_const_for_grazing,
                            preference_for_nanophytoplankton,
                            preference_for_diatoms,
                            preference_for_POC,
                            preference_for_microzooplankton,
                            food_threshold_for_zooplankton,
                            specific_food_thresholds_for_microzooplankton,
                            specific_food_thresholds_for_mesozooplankton,
                            zooplankton_quadratic_mortality,
                            zooplankton_linear_mortality,
                            half_saturation_const_for_mortality,
                            fraction_of_calcite_not_dissolving_in_guts,
                            FeC_ratio_of_zooplankton,
                            FeZ_redfield_ratio, 


                            remineralisation_rate_of_DOC,
                            half_saturation_const_for_DOC_remin,
                            NO3_half_saturation_const_for_DOC_remin,
                            NH4_half_saturation_const_for_DOC_remin,
                            PO4_half_saturation_const_for_DOC_remin,
                            Fe_half_saturation_const_for_DOC_remin,
                            aggregation_rate_of_DOC_to_POC_1,
                            aggregation_rate_of_DOC_to_POC_2,
                            aggregation_rate_of_DOC_to_GOC_3,
                            aggregation_rate_of_DOC_to_POC_4,
                            aggregation_rate_of_DOC_to_POC_5,


                            degradation_rate_of_POC,
                            sinking_speed_of_POC,
                            min_sinking_speed_of_GOC,
                            sinking_speed_of_dust,
                            aggregation_rate_of_POC_to_GOC_6,
                            aggregation_rate_of_POC_to_GOC_7,
                            aggregation_rate_of_POC_to_GOC_8,
                            aggregation_rate_of_POC_to_GOC_9,
                            min_scavenging_rate_of_iron,
                            slope_of_scavenging_rate_of_iron,
                            scavenging_rate_of_iron_by_dust,
                            dissolution_rate_of_calcite,
                            exponent_in_the_dissolution_rate_of_calcite,
                            proportion_of_the_most_labile_phase_in_PSi,
                            slow_dissolution_rate_of_PSi,
                            fast_dissolution_rate_of_PSi,


                            max_nitrification_rate,
                            half_sat_const_for_denitrification1,
                            half_sat_const_for_denitrification2,
                            total_concentration_of_iron_ligands,
                            max_rate_of_nitrogen_fixation,
                            Fe_half_saturation_constant_of_nitrogen_fixation,
                            photosynthetic_parameter_of_nitrogen_fixation,
                            iron_concentration_in_sea_ice,
                            max_sediment_flux_of_Fe,
                            solubility_of_iron_in_dust,
                            OC_for_ammonium_based_processes,
                            OC_ratio_of_nitrification,
                            CN_ratio_of_ammonification,
                            CN_ratio_of_denitrification,
                            NC_redfield_ratio,
                            PC_redfield_ratio,
                            rain_ratio_parameter,
                            bacterial_reference,

                            NC_stoichiometric_ratio_of_dentitrification,
                            NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER,
                            dissolution_rate_of_silicon,
                            coefficient_of_bacterial_uptake_of_iron_in_POC,
                            coefficient_of_bacterial_uptake_of_iron_in_GOC,
                            max_FeC_ratio_of_bacteria,
                            Fe_half_saturation_const_for_PLACEHOLDER,    #not sure what this should be called
                            proportion_of_sinking_grazed_shells,

                            mixed_layer_depth,
                            euphotic_layer_depth,
                            yearly_maximum_silicate,
                            dust_deposition,
                            vertical_diffusivity,
                            carbonate_sat_ratio,

                          sinking_velocities)
    end
    
end

"""
    PISCES(; grid, # finally the function
                   # now you can finally put the values here
                   growth_rate_at_zero :: FT = 0.6 / day,                                 # 1/d,
                   growth_rate_reference_for_light_limitation :: FT = 1.0/ day,           # 1/d
                   basal_respiration_rate :: FT = 0.033 / day,                             # 1/d
                   temperature_sensitivity_of_growth :: FT = 1.066,
                   initial_slope_of_PI_curve :: PD = (P = 2/day, D = 2/day),        #(Wm⁻²)⁻¹d⁻¹  
                   exudation_of_DOC :: PD = (P = 0.05, D = 0.05),  
                   absorption_in_the_blue_part_of_light :: PD = (P = 2.1, D = 1.6),
                   absorption_in_the_green_part_of_light :: PD = (P = 0.42, D = 0.69),
                   absorption_in_the_red_part_of_light :: PD = (P = 0.4, D = 0.7),
                   min_half_saturation_const_for_phosphate :: PD = (P = 0.8, D = 2.4),     #nmolPL⁻¹    
                   min_half_saturation_const_for_ammonium :: PD = (P = 0.013, D = 0.039),  #μmolNL⁻¹
                   min_half_saturation_const_for_nitrate :: PD = (P = 0.13, D =0.39),     #μmolNL⁻¹
                   min_half_saturation_const_for_silicate :: FT = 1.0,            #μmolSiL⁻¹
                   parameter_for_half_saturation_const :: FT = 16.6,            #μmolSiL⁻¹
                   parameter_for_SiC :: OT = (one = 2.0, two = 20.0),                           #μmolSiL⁻¹
                   min_half_saturation_const_for_iron_uptake :: PD = (P = 1.0, D = 3.0),   #nmolFeL⁻¹
                   size_ratio_of_phytoplankton :: PD = (P = 3.0, D = 3.0),
                   optimal_SiC_uptake_ratio_of_diatoms :: FT = 0.159,       #molSi/(mol C)
                   optimal_iron_quota :: PD = (P = 7.0, D = 7.0),               #μmolFe/(mol C)
                   max_iron_quota :: PD = (P = 40.0, D = 40.0),                  #μmolFe/(mol C)
                   phytoplankton_mortality_rate :: PD = (P = 0.01/day, D = 0.01/day),
                   min_quadratic_mortality_of_phytoplankton :: FT = 0.01 / day,   #1/(d mol C)
                   max_quadratic_mortality_of_diatoms :: FT = 0.03 / day,         #1/(d mol C)
                   max_ChlC_ratios_of_phytoplankton :: PD = (P = 0.033, D = 0.05),  #mg Chl/(mg C)
                   min_ChlC_ratios_of_phytoplankton :: FT = 0.0033,    #mg Chl/(mg C)
                   threshold_concentration_for_size_dependency :: PD = (P = 1.0, D = 1.0),  #μmolCL⁻¹
                   mean_residence_time_of_phytoplankton_in_unlit_mixed_layer :: PD = (P = 3.0/day, D = 4.0/day), #/day
    
                   latitude :: FT = -1.0, #still to be changed - this is temporary 
                   length_of_day :: FT = 1.0, #temporary parameter for day length

                   temperature_sensitivity_term :: ZM = (Z = 1.079, M = 1.079),   
                   max_growth_efficiency_of_zooplankton :: ZM = (Z = 0.3, M = 0.35),
                   non_assimilated_fraction :: ZM = (Z = 0.3, M = 0.3),
                   excretion_as_DOM :: ZM = (Z = 0.6, M = 0.6),
                   max_grazing_rate :: ZM = (Z = 3.0/day, M = 0.75/day),                       #1/d
                   flux_feeding_rate :: FT = 2.0e-3,                                #(m mol L⁻¹)⁻¹
                   half_saturation_const_for_grazing :: ZM = (Z = 20.0, M = 20.0),               #μmolCL⁻¹
                   preference_for_nanophytoplankton :: ZM = (Z = 1.0, M = 0.3),
                   preference_for_diatoms :: ZM = (Z = 0.5, M = 1.0),
                   preference_for_POC :: ZM= (Z = 0.1, M = 0.3),
                   preference_for_microzooplankton :: FT = 1.0,
                   food_threshold_for_zooplankton :: ZM = (Z = 0.3, M = 0.3),                  #μmolCL⁻¹
                   specific_food_thresholds_for_microzooplankton :: FT = 0.001,        #μmolCL⁻¹
                   specific_food_thresholds_for_mesozooplankton :: FT = 0.001,         #μmolCL⁻¹
                   zooplankton_quadratic_mortality :: ZM = (Z = 0.004/day, M = 0.03/day),       #(μmolCL⁻¹)⁻¹d⁻¹
                   zooplankton_linear_mortality :: ZM = (Z = 0.03/day, M = 0.005/day),           #1/d
                   half_saturation_const_for_mortality :: FT = 0.2,                     #μmolCL⁻¹
                   fraction_of_calcite_not_dissolving_in_guts :: ZM = (Z = 0.5, M = 0.75),
                   FeC_ratio_of_zooplankton :: FT = 10.0,                                  #μmolFe molC⁻¹
                   FeZ_redfield_ratio :: FT = 3.0,             #μmolFe molC⁻¹
   
   
                   remineralisation_rate_of_DOC :: FT = 0.3 / day,                 #1/d
                   half_saturation_const_for_DOC_remin :: FT = 417.0,                #μmolCL⁻¹
                   NO3_half_saturation_const_for_DOC_remin :: FT = 0.03,            #μmolNL⁻¹
                   NH4_half_saturation_const_for_DOC_remin :: FT = 0.003,           #μmolNL⁻¹
                   PO4_half_saturation_const_for_DOC_remin :: FT = 0.003,       #μmolPL⁻¹
                   Fe_half_saturation_const_for_DOC_remin :: FT = 0.01,         #μmolFeL⁻¹
                   aggregation_rate_of_DOC_to_POC_1 :: FT = 0.37 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_2 :: FT = 102.0 / day,           #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_GOC_3 :: FT = 3530.0 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_4 :: FT = 5095.0 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_5 :: FT = 114.0 / day,           #(μmolCL⁻¹)⁻¹d⁻¹
   
   
                   degradation_rate_of_POC :: FT = 0.025 / day,             #1/d
                   sinking_speed_of_POC :: FT = 2.0 / day,                    #md⁻¹
                   min_sinking_speed_of_GOC :: FT = 30.0 / day,               #md⁻¹
                   sinking_speed_of_dust :: FT = 2.0,                         #ms⁻¹
                   aggregation_rate_of_POC_to_GOC_6 :: FT = 25.9 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_7 :: FT = 4452 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_8 :: FT = 3.3 / day,      #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_9 :: FT = 47.1 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   min_scavenging_rate_of_iron :: FT = 3.0e-5 / day,          #1/d
                   slope_of_scavenging_rate_of_iron :: FT = 0.005 / day,    #d⁻¹μmol⁻¹L
                   scavenging_rate_of_iron_by_dust :: FT = 150.0 / day,       #d⁻¹mg⁻¹L
                   dissolution_rate_of_calcite :: FT = 0.197 / day,         #1/d
                   exponent_in_the_dissolution_rate_of_calcite :: FT = 1.0,
                   proportion_of_the_most_labile_phase_in_PSi :: FT = 0.5,
                   slow_dissolution_rate_of_PSi :: FT = 0.003 / day,        #1/d
                   fast_dissolution_rate_of_PSi :: FT = 0.025 / day,        #1/d
   
   
                   max_nitrification_rate :: FT = 0.05 / day,                           #1/d
                   half_sat_const_for_denitrification1 :: FT = 1.0,                       #μmolO₂L⁻¹
                   half_sat_const_for_denitrification2 :: FT = 6.0,                       #μmolO₂L⁻¹
                   total_concentration_of_iron_ligands :: FT = 0.6,                     #nmolL⁻¹
                   max_rate_of_nitrogen_fixation :: FT = 0.013,                         #μmolNL⁻¹d⁻¹
                   Fe_half_saturation_constant_of_nitrogen_fixation :: FT = 0.1,        #nmolFeL⁻¹
                   photosynthetic_parameter_of_nitrogen_fixation :: FT = 50.0,            #Wm⁻²
                   iron_concentration_in_sea_ice :: FT = 15.0,                            #nmolFeL⁻¹   
                   max_sediment_flux_of_Fe :: FT = 2.0 / day,                             #μmolFem⁻²d⁻¹
                   solubility_of_iron_in_dust :: FT = 0.02,
                   OC_for_ammonium_based_processes :: FT = 133/122,                     #molO₂(mol C)⁻¹
                   OC_ratio_of_nitrification :: FT = 32/122,                            #molO₂(mol C)⁻¹
                   CN_ratio_of_ammonification :: FT = 3/5,                              #molN(mol C)⁻¹
                   CN_ratio_of_denitrification :: FT = 105/16,                          #molN(mol C)⁻¹
                   NC_redfield_ratio :: FT = 16/122,    
                   PC_redfield_ratio :: = 1/122,                                #molN(mol C)⁻¹
                   rain_ratio_parameter :: FT = 0.3,
                   bacterial_reference :: FT = 1.0,     #Not sure if this is what its called : denoted Bact_ref in paper

                   NC_stoichiometric_ratio_of_dentitrification :: FT = 0.86,
                   NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER :: FT = 0.0,     #again not sure what this is called
                   dissolution_rate_of_silicon :: FT = 1.0,
                   coefficient_of_bacterial_uptake_of_iron_in_POC :: FT = 0.5,
                   coefficient_of_bacterial_uptake_of_iron_in_GOC :: FT = 0.5,
                   max_FeC_ratio_of_bacteria :: FT = 10.0e-6,     #or 6
                   Fe_half_saturation_const_for_PLACEHOLDER :: FT = 2.5e-10, #or 2.5e-10    #not sure what this should be called
                   proportion_of_sinking_grazed_shells :: ZM = (Z = 0.3, M = 0.3),  # 0.3 for both? not sure

                   mixed_layer_depth :: CF = ConstantField(100),
                   euphotic_layer_depth :: CF = ConstantField(50),
                   vertical_diffusivity :: CF  = ConstantField(1),
                   yearly_maximum_silicate :: FT = 1.0,
                   dust_deposition :: FT = 1.0,

                  surface_photosynthetically_active_radiation = default_surface_PAR,

                  light_attenuation_model::LA =
                    TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                            surface_PAR = surface_photosynthetically_active_radiation),

                  # just keep all this stuff for now but you can ignore it
                  sediment_model::S = nothing,

                  sinking_speeds = (POC = 0.0, GOC = 0.0, SFe = 0.0, BFe = 0.0, PSi = 0.0, CaCO₃ = 0.0),  #change all 1.0s to w_GOC
                  
                  carbonate_sat_ratio :: ZF = ZeroField(),
                  open_bottom::Bool = true,

                  scale_negatives = false,

                  particles::P = nothing,
                  modifiers::M = nothing)

Construct an instance of the [PISCES](@ref PISCES) biogeochemical model. 

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to calculate sinking
- `parameter_1`...: PISCES parameter values
- `surface_photosynthetically_active_radiation`: funciton (or array in the future) for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `sinking_speed`: named tuple of constant sinking, of fields (i.e. `ZFaceField(...)`) for any tracers which sink (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
- `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME

julia> using Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> model = PISCES(; grid)
PISCES{Float64} ... # we can fix this later
```
"""
function PISCES(; grid, # finally the function
                   # now you can finally put the values here
                   growth_rate_at_zero :: FT = 0.6 / day,                                 # 1/d,
                   growth_rate_reference_for_light_limitation :: FT = 1.0/ day,           # 1/d
                   basal_respiration_rate :: FT = 0.033 / day,                             # 1/d
                   temperature_sensitivity_of_growth :: FT = 1.066,
                   initial_slope_of_PI_curve :: PD = (P = 2/day, D = 2/day),        #(Wm⁻²)⁻¹d⁻¹  
                   exudation_of_DOC :: PD = (P = 0.05, D = 0.05),  
                   absorption_in_the_blue_part_of_light :: PD = (P = 2.1, D = 1.6),
                   absorption_in_the_green_part_of_light :: PD = (P = 0.42, D = 0.69),
                   absorption_in_the_red_part_of_light :: PD = (P = 0.4, D = 0.7),
                   min_half_saturation_const_for_phosphate :: PD = (P = 0.8, D = 2.4),     #nmolPL⁻¹    
                   min_half_saturation_const_for_ammonium :: PD = (P = 0.013, D = 0.039),  #μmolNL⁻¹
                   min_half_saturation_const_for_nitrate :: PD = (P = 0.13, D =0.39),     #μmolNL⁻¹
                   min_half_saturation_const_for_silicate :: FT = 1.0,            #μmolSiL⁻¹
                   parameter_for_half_saturation_const :: FT = 16.6,            #μmolSiL⁻¹
                   parameter_for_SiC :: OT = (one = 2.0, two = 20.0),                           #μmolSiL⁻¹
                   min_half_saturation_const_for_iron_uptake :: PD = (P = 1.0, D = 3.0),   #nmolFeL⁻¹
                   size_ratio_of_phytoplankton :: PD = (P = 3.0, D = 3.0),
                   optimal_SiC_uptake_ratio_of_diatoms :: FT = 0.159,       #molSi/(mol C)
                   optimal_iron_quota :: PD = (P = 7.0, D = 7.0),               #mmolFe/(mol C)
                   max_iron_quota :: PD = (P = 40.0, D = 40.0),                  #molFe/(mol C)
                   phytoplankton_mortality_rate :: PD = (P = 0.01/day, D = 0.01/day),
                   min_quadratic_mortality_of_phytoplankton :: FT = 0.01 / day,   #1/(d mol C)
                   max_quadratic_mortality_of_diatoms :: FT = 0.03 / day,         #1/(d mol C)
                   max_ChlC_ratios_of_phytoplankton :: PD = (P = 0.033, D = 0.05),  #mg Chl/(mg C)
                   min_ChlC_ratios_of_phytoplankton :: FT = 0.0033,    #mg Chl/(mg C)
                   threshold_concentration_for_size_dependency :: PD = (P = 1.0, D = 1.0),  #μmolCL⁻¹
                   mean_residence_time_of_phytoplankton_in_unlit_mixed_layer :: PD = (P = 3days, D = 4days), #day
    
                   latitude :: FT = -1.0, #still to be changed - this is temporary 
                   length_of_day :: FT = 1.0, #temporary parameter for day length

                   temperature_sensitivity_term :: ZM = (Z = 1.079, M = 1.079),   
                   max_growth_efficiency_of_zooplankton :: ZM = (Z = 0.3, M = 0.35),
                   non_assimilated_fraction :: ZM = (Z = 0.3, M = 0.3),
                   excretion_as_DOM :: ZM = (Z = 0.6, M = 0.6),
                   max_grazing_rate :: ZM = (Z = 3.0/day, M = 0.75/day),                       #1/d
                   flux_feeding_rate :: FT = 2.0e-3,                                #(m mol L⁻¹)⁻¹
                   half_saturation_const_for_grazing :: ZM = (Z = 20.0, M = 20.0),               #μmolCL⁻¹
                   preference_for_nanophytoplankton :: ZM = (Z = 1.0, M = 0.3),
                   preference_for_diatoms :: ZM = (Z = 0.5, M = 1.0),
                   preference_for_POC :: ZM= (Z = 0.1, M = 0.3),
                   preference_for_microzooplankton :: FT = 1.0,
                   food_threshold_for_zooplankton :: ZM = (Z = 0.3, M = 0.3),                  #μmolCL⁻¹
                   specific_food_thresholds_for_microzooplankton :: FT = 0.001,        #μmolCL⁻¹
                   specific_food_thresholds_for_mesozooplankton :: FT = 0.001,         #μmolCL⁻¹
                   zooplankton_quadratic_mortality :: ZM = (Z = 0.004/day, M = 0.03/day),       #(μmolCL⁻¹)⁻¹d⁻¹
                   zooplankton_linear_mortality :: ZM = (Z = 0.03/day, M = 0.005/day),           #1/d
                   half_saturation_const_for_mortality :: FT = 0.2,                     #μmolCL⁻¹
                   fraction_of_calcite_not_dissolving_in_guts :: ZM = (Z = 0.5, M = 0.75),
                   FeC_ratio_of_zooplankton :: FT = 10.0e-3,                                  #mmolFe molC⁻¹
                   FeZ_redfield_ratio :: FT = 3.0e-3,             #mmolFe molC⁻¹, remove this, is actually FeC_ratio_of_zooplankton
   
   
                   remineralisation_rate_of_DOC :: FT = 0.3 / day,                 #1/d
                   half_saturation_const_for_DOC_remin :: FT = 417.0,                #μmolCL⁻¹
                   NO3_half_saturation_const_for_DOC_remin :: FT = 0.03,            #μmolNL⁻¹
                   NH4_half_saturation_const_for_DOC_remin :: FT = 0.003,           #μmolNL⁻¹
                   PO4_half_saturation_const_for_DOC_remin :: FT = 0.003,       #μmolPL⁻¹
                   Fe_half_saturation_const_for_DOC_remin :: FT = 0.01,         #μmolFeL⁻¹
                   aggregation_rate_of_DOC_to_POC_1 :: FT = 0.37 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_2 :: FT = 102.0 / day,           #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_GOC_3 :: FT = 3530.0 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_4 :: FT = 5095.0 / day,          #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_DOC_to_POC_5 :: FT = 114.0 / day,           #(μmolCL⁻¹)⁻¹d⁻¹
   
   
                   degradation_rate_of_POC :: FT = 0.025 / day,             #1/d
                   sinking_speed_of_POC :: FT = 2.0 / day,                    #md⁻¹
                   min_sinking_speed_of_GOC :: FT = 30.0 / day,               #md⁻¹
                   sinking_speed_of_dust :: FT = 2.0,                         #ms⁻¹
                   aggregation_rate_of_POC_to_GOC_6 :: FT = 25.9 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_7 :: FT = 4452.0 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_8 :: FT = 3.3 / day,      #(μmolCL⁻¹)⁻¹d⁻¹
                   aggregation_rate_of_POC_to_GOC_9 :: FT = 47.1 / day,     #(μmolCL⁻¹)⁻¹d⁻¹
                   min_scavenging_rate_of_iron :: FT = 3.0e-5 / day,          #1/d
                   slope_of_scavenging_rate_of_iron :: FT = 0.005 / day,    #d⁻¹μmol⁻¹L
                   scavenging_rate_of_iron_by_dust :: FT = 150.0 / day,       #d⁻¹mg⁻¹L
                   dissolution_rate_of_calcite :: FT = 0.197 / day,         #1/d
                   exponent_in_the_dissolution_rate_of_calcite :: FT = 1.0,
                   proportion_of_the_most_labile_phase_in_PSi :: FT = 0.5,
                   slow_dissolution_rate_of_PSi :: FT = 0.003 / day,        #1/d
                   fast_dissolution_rate_of_PSi :: FT = 0.025 / day,        #1/d
   
   
                   max_nitrification_rate :: FT = 0.05 / day,                           #1/d
                   half_sat_const_for_denitrification1 :: FT = 1.0,                       #μmolO₂L⁻¹
                   half_sat_const_for_denitrification2 :: FT = 6.0,                       #μmolO₂L⁻¹
                   total_concentration_of_iron_ligands :: FT = 0.6,                     #nmolL⁻¹
                   max_rate_of_nitrogen_fixation :: FT = 0.013,                         #μmolNL⁻¹d⁻¹
                   Fe_half_saturation_constant_of_nitrogen_fixation :: FT = 0.1,        #nmolFeL⁻¹
                   photosynthetic_parameter_of_nitrogen_fixation :: FT = 50.0,            #Wm⁻²
                   iron_concentration_in_sea_ice :: FT = 15.0,                            #nmolFeL⁻¹   
                   max_sediment_flux_of_Fe :: FT = 2.0 / day,                             #μmolFem⁻²d⁻¹
                   solubility_of_iron_in_dust :: FT = 0.02,
                   OC_for_ammonium_based_processes :: FT = 133/122,                     #molO₂(mol C)⁻¹
                   OC_ratio_of_nitrification :: FT = 32/122,                            #molO₂(mol C)⁻¹
                   CN_ratio_of_ammonification :: FT = 3/5,                              #molN(mol C)⁻¹
                   CN_ratio_of_denitrification :: FT = 105/16,                          #molN(mol C)⁻¹
                   NC_redfield_ratio :: FT = 16/122, 
                   PC_redfield_ratio :: FT = 1/122,                                   #molN(mol C)⁻¹
                   rain_ratio_parameter :: FT = 0.3,
                   bacterial_reference :: FT = 1.0,     #Not sure if this is what its called : denoted Bact_ref in paper

                   NC_stoichiometric_ratio_of_dentitrification :: FT = 0.86,
                   NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER :: FT = 0.0,     #again not sure what this is called
                   dissolution_rate_of_silicon :: FT = 1.0,
                   coefficient_of_bacterial_uptake_of_iron_in_POC :: FT = 0.5,
                   coefficient_of_bacterial_uptake_of_iron_in_GOC :: FT = 0.5,
                   max_FeC_ratio_of_bacteria :: FT = 10.0e-6,     #or 6
                   Fe_half_saturation_const_for_PLACEHOLDER :: FT = 2.5e-10, #or 2.5e-10    #not sure what this should be called
                   proportion_of_sinking_grazed_shells :: ZM = (Z = 0.3, M = 0.3),  # 0.3 for both? not sure

                   mixed_layer_depth :: CF = ConstantField(100),
                   euphotic_layer_depth :: CF = ConstantField(50),
                   vertical_diffusivity :: CF  = ConstantField(1),
                   yearly_maximum_silicate :: FT = 1.0,
                   dust_deposition :: FT = 1.0,

                  surface_photosynthetically_active_radiation = default_surface_PAR,

                  light_attenuation_model::LA =
                    TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                            surface_PAR = surface_photosynthetically_active_radiation),

                  # just keep all this stuff for now but you can ignore it
                  sediment_model::S = nothing,

                  sinking_speeds = (POC = 0.0, GOC = 0.0, SFe = 0.0, BFe = 0.0, PSi = 0.0, CaCO₃ = 0.0),  #change all 1.0s to w_GOC
                  
                  carbonate_sat_ratio :: ZF = ZeroField(),
                  open_bottom::Bool = true,

                  scale_negatives = false,

                  particles::P = nothing,
                  modifiers::M = nothing) where {FT, PD, ZM, OT, LA, S, P, M, CF, ZF}

    if !isnothing(sediment_model) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    underlying_biogeochemistry = PISCES(growth_rate_at_zero,
                                        growth_rate_reference_for_light_limitation,
                                        basal_respiration_rate,
                                        temperature_sensitivity_of_growth,
                                        initial_slope_of_PI_curve,
                                        exudation_of_DOC,
                                        absorption_in_the_blue_part_of_light,
                                        absorption_in_the_green_part_of_light,
                                        absorption_in_the_red_part_of_light,
                                        min_half_saturation_const_for_phosphate,
                                        min_half_saturation_const_for_ammonium,
                                        min_half_saturation_const_for_nitrate,
                                        min_half_saturation_const_for_silicate,
                                        parameter_for_half_saturation_const,
                                        parameter_for_SiC,
                                        min_half_saturation_const_for_iron_uptake,
                                        size_ratio_of_phytoplankton,
                                        optimal_SiC_uptake_ratio_of_diatoms,
                                        optimal_iron_quota,
                                        max_iron_quota,
                                        phytoplankton_mortality_rate,
                                        min_quadratic_mortality_of_phytoplankton,
                                        max_quadratic_mortality_of_diatoms,
                                        max_ChlC_ratios_of_phytoplankton,
                                        min_ChlC_ratios_of_phytoplankton,
                                        threshold_concentration_for_size_dependency,
                                        mean_residence_time_of_phytoplankton_in_unlit_mixed_layer,

                                        latitude,
                                        length_of_day,
                                        
                                        
                                        temperature_sensitivity_term,
                                        max_growth_efficiency_of_zooplankton,
                                        non_assimilated_fraction,
                                        excretion_as_DOM,
                                        max_grazing_rate,
                                        flux_feeding_rate,
                                        half_saturation_const_for_grazing,
                                        preference_for_nanophytoplankton,
                                        preference_for_diatoms,
                                        preference_for_POC,
                                        preference_for_microzooplankton,
                                        food_threshold_for_zooplankton,
                                        specific_food_thresholds_for_microzooplankton,
                                        specific_food_thresholds_for_mesozooplankton,
                                        zooplankton_quadratic_mortality,
                                        zooplankton_linear_mortality,
                                        half_saturation_const_for_mortality,
                                        fraction_of_calcite_not_dissolving_in_guts,
                                        FeC_ratio_of_zooplankton,
                                        FeZ_redfield_ratio,


                                        remineralisation_rate_of_DOC,
                                        half_saturation_const_for_DOC_remin,
                                        NO3_half_saturation_const_for_DOC_remin,
                                        NH4_half_saturation_const_for_DOC_remin,
                                        PO4_half_saturation_const_for_DOC_remin,
                                        Fe_half_saturation_const_for_DOC_remin,
                                        aggregation_rate_of_DOC_to_POC_1,
                                        aggregation_rate_of_DOC_to_POC_2,
                                        aggregation_rate_of_DOC_to_GOC_3,
                                        aggregation_rate_of_DOC_to_POC_4,
                                        aggregation_rate_of_DOC_to_POC_5,


                                        degradation_rate_of_POC,
                                        sinking_speed_of_POC,
                                        min_sinking_speed_of_GOC,
                                        sinking_speed_of_dust,
                                        aggregation_rate_of_POC_to_GOC_6,
                                        aggregation_rate_of_POC_to_GOC_7,
                                        aggregation_rate_of_POC_to_GOC_8,
                                        aggregation_rate_of_POC_to_GOC_9,
                                        min_scavenging_rate_of_iron,
                                        slope_of_scavenging_rate_of_iron,
                                        scavenging_rate_of_iron_by_dust,
                                        dissolution_rate_of_calcite,
                                        exponent_in_the_dissolution_rate_of_calcite,
                                        proportion_of_the_most_labile_phase_in_PSi,
                                        slow_dissolution_rate_of_PSi,
                                        fast_dissolution_rate_of_PSi,


                                        max_nitrification_rate,
                                        half_sat_const_for_denitrification1,
                                        half_sat_const_for_denitrification2,
                                        total_concentration_of_iron_ligands,
                                        max_rate_of_nitrogen_fixation,
                                        Fe_half_saturation_constant_of_nitrogen_fixation,
                                        photosynthetic_parameter_of_nitrogen_fixation,
                                        iron_concentration_in_sea_ice,
                                        max_sediment_flux_of_Fe,
                                        solubility_of_iron_in_dust,
                                        OC_for_ammonium_based_processes,
                                        OC_ratio_of_nitrification,
                                        CN_ratio_of_ammonification,
                                        CN_ratio_of_denitrification,
                                        NC_redfield_ratio,
                                        PC_redfield_ratio,
                                        rain_ratio_parameter,
                                        bacterial_reference,

                                        NC_stoichiometric_ratio_of_dentitrification,
                                        NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER,
                                        dissolution_rate_of_silicon,
                                        coefficient_of_bacterial_uptake_of_iron_in_POC,
                                        coefficient_of_bacterial_uptake_of_iron_in_GOC,
                                        max_FeC_ratio_of_bacteria,
                                        Fe_half_saturation_const_for_PLACEHOLDER,    #not sure what this should be called
                                        proportion_of_sinking_grazed_shells,

                                        mixed_layer_depth,
                                        euphotic_layer_depth,
                                        yearly_maximum_silicate,
                                        dust_deposition,


                                        vertical_diffusivity,
                                        carbonate_sat_ratio,

                                        sinking_velocities)

    if scale_negatives
        scaler = ScaleNegativeTracers(underlying_biogeochemistry, grid; invalid_fill_value)
        if isnothing(modifiers)
            modifiers = scaler
        elseif modifiers isa Tuple
            modifiers = (modifiers..., scaler)
        else
            modifiers = (modifiers, scaler)
        end
    end

    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation = light_attenuation_model, 
                           sediment = sediment_model, 
                           particles,
                           modifiers)
end

@inline biogeochemical_auxiliary_fields(bgc::PISCES) = (zₘₓₗ = bgc.mixed_layer_depth, zₑᵤ = bgc.euphotic_layer_depth, Si̅ = bgc.yearly_maximum_silicate, D_dust = bgc.dust_deposition, Ω = bgc.carbonate_sat_ratio)

@inline required_biogeochemical_tracers(::PISCES) = (:P, :D, :Z, :M, :Pᶜʰˡ, :Dᶜʰˡ, :Pᶠᵉ, :Dᶠᵉ, :Dˢⁱ, :DOC, :POC, :GOC, :SFe, :BFe, :PSi, :NO₃, :NH₄, :PO₄, :Fe, :Si, :CaCO₃, :DIC, :Alk, :O₂, :T) # list all the parameters here, also if you need T and S put them here too

@inline required_biogeochemical_auxiliary_fields(::PISCES) = (:zₘₓₗ, :zₑᵤ, :Si̅, :D_dust, :Ω, :PAR, :PAR¹, :PAR², :PAR³, )

# for sinking things like POM this is how we tell oceananigans ther sinking speed
@inline function biogeochemical_drift_velocity(bgc::PISCES, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end



# don't worry about this for now
adapt_structure(to, pisces::PISCES) =
    PISCES(adapt(to, pisces.parameter_1),
           adapt(to, pisces.sinking_velocities))

# you can updatye these if you want it to have a pretty way of showing uyou its a pisces model
summary(::PISCES{FT}) where {FT} = string("PISCES{$FT}") 

show(io::IO, model::PISCES) where {FT, B, W, PD, ZM, OT}  = print(io, string("Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model")) # maybe add some more info here

@inline maximum_sinking_velocity(bgc::PISCES) = maximum(abs, bgc.sinking_velocities.bPOM.w) # might need ot update this for wghatever the fastest sinking pareticles are

# write most of the code here (i.e. make a file falled phytoplankton.jl and then include it here)
include("phytoplankton.jl")
include("calcite.jl")
include("carbonate_system.jl")
include("DOC.jl")
include("iron_in_particles.jl")
include("iron.jl")
include("nitrates_ammonium.jl")
include("oxygen.jl")
include("phosphates.jl")
include("POC_and_GOC.jl")
include("psi.jl")
include("si.jl")
include("zooplankton.jl")

# to work with the sediment model we need to tell in the redfield ratio etc. of some things, but for now we can ignore
@inline redfield(i, j, k, val_tracer_name, bgc::PISCES, tracers) = NaN

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline remineralisation_receiver(::PISCES) = :NH₄

# this is for positivity preservation, if you can work it out it would be great, I don't think PISCES conserves C but probably does Nitrogen
@inline conserved_tracers(::PISCES) = NaN

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here
end # module
