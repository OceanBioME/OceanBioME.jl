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


using Oceananigans: KernelFunctionOperation
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField, ConstantField, Center, Face

using OceanBioME.Light: MultiBandPhotosyntheticallyActiveRadiation, default_surface_PAR, compute_euphotic_depth!
using OceanBioME: setup_velocity_fields, show_sinking_velocities, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.BoxModels: BoxModel

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Grids: φnodes, RectilinearGrid

import OceanBioME: redfield, conserved_tracers, maximum_sinking_velocity, chlorophyll


import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity,
                                     biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt
import Base: show, summary

import OceanBioME.Models.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

include("common.jl")

struct PISCES{FT, PD, ZM, OT, W, CF, ZF, LA, FFMLD, FFEU} <: AbstractContinuousFormBiogeochemistry
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

    latitude :: LA

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
    Fe_half_saturation_const_for_Bacteria :: FT    #not sure what this should be called

    mixed_layer_shear :: FT
    background_shear :: FT

    mixed_layer_depth :: FFMLD
    euphotic_depth :: FFEU
    yearly_maximum_silicate :: CF
    dust_deposition :: ZF

    vertical_diffusivity :: CF 
    carbonate_sat_ratio :: ZF

    sinking_velocities :: W
end

"""
    PISCES(; grid,
 

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
PISCES{Float64} 
 Light attenuation: Two-band light attenuation model (Float64)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing
```
"""
function PISCES(; grid,
                   growth_rate_at_zero = 0.6 / day,                                 # 1/second
                   growth_rate_reference_for_light_limitation = 1.0/ day,           # 1/second
                   basal_respiration_rate = 0.033 / day,                             # 1/second
                   temperature_sensitivity_of_growth = 1.066,
                   initial_slope_of_PI_curve = (P = 2/day, D = 2/day),        #(Wm⁻²)⁻¹s⁻¹  
                   exudation_of_DOC = (P = 0.05, D = 0.05),  
                   absorption_in_the_blue_part_of_light = (P = 2.1, D = 1.6),
                   absorption_in_the_green_part_of_light = (P = 0.42, D = 0.69),
                   absorption_in_the_red_part_of_light = (P = 0.4, D = 0.7),
                   min_half_saturation_const_for_phosphate = (P = 0.8, D = 2.4),     #nmolPL⁻¹    
                   min_half_saturation_const_for_ammonium = (P = 0.013, D = 0.039),  #μmolNL⁻¹
                   min_half_saturation_const_for_nitrate = (P = 0.13, D =0.39),     #μmolNL⁻¹
                   min_half_saturation_const_for_silicate = 1.0,            #μmolSiL⁻¹
                   parameter_for_half_saturation_const = 16.6,            #μmolSiL⁻¹
                   parameter_for_SiC = (one = 2.0, two = 20.0),                           #μmolSiL⁻¹
                   min_half_saturation_const_for_iron_uptake = (P = 1.0, D = 3.0),   #nmolFeL⁻¹
                   size_ratio_of_phytoplankton = (P = 3.0, D = 3.0),
                   optimal_SiC_uptake_ratio_of_diatoms = 0.159,       #molSi/(mol C)
                   optimal_iron_quota = (P = 7.0e-3, D = 7.0e-3),               #mmolFe/(mol C)
                   max_iron_quota = (P = 40.0e-3, D = 40.0e-3),                  #molFe/(mol C)
                   phytoplankton_mortality_rate = (P = 0.01/day, D = 0.01/day), #1/second
                   min_quadratic_mortality_of_phytoplankton = 0.01 / day,   #1/(d mol C)
                   max_quadratic_mortality_of_diatoms = 0.03 / day,         #1/(d mol C)
                   max_ChlC_ratios_of_phytoplankton = (P = 0.033, D = 0.05),  #mg Chl/(mg C)
                   min_ChlC_ratios_of_phytoplankton = 0.0033,    #mg Chl/(mg C)
                   threshold_concentration_for_size_dependency = (P = 1.0, D = 1.0),  #μmolCL⁻¹
                   mean_residence_time_of_phytoplankton_in_unlit_mixed_layer = (P = 3days, D = 4days), #seconds
    
                   latitude = 45,

                   temperature_sensitivity_term = (Z = 1.079, M = 1.079),   
                   max_growth_efficiency_of_zooplankton = (Z = 0.3, M = 0.35),
                   non_assimilated_fraction = (Z = 0.3, M = 0.3),
                   excretion_as_DOM = (Z = 0.6, M = 0.6),
                   max_grazing_rate = (Z = 3.0/day, M = 0.75/day),                       #1/second
                   flux_feeding_rate = 2.0e-3,                                #(m mol L⁻¹)⁻¹
                   half_saturation_const_for_grazing = (Z = 20.0, M = 20.0),               #μmolCL⁻¹
                   preference_for_nanophytoplankton = (Z = 1.0, M = 0.3),
                   preference_for_diatoms = (Z = 0.5, M = 1.0),
                   preference_for_POC= (Z = 0.1, M = 0.3),
                   preference_for_microzooplankton = 1.0,
                   food_threshold_for_zooplankton = (Z = 0.3, M = 0.3),                  #μmolCL⁻¹
                   specific_food_thresholds_for_microzooplankton = 0.001,        #μmolCL⁻¹
                   specific_food_thresholds_for_mesozooplankton = 0.001,         #μmolCL⁻¹
                   zooplankton_quadratic_mortality = (Z = 0.004/day, M = 0.03/day),       #(μmolCL⁻¹)⁻¹s⁻¹
                   zooplankton_linear_mortality = (Z = 0.03/day, M = 0.005/day),           #1/second
                   half_saturation_const_for_mortality = 0.2,                     #μmolCL⁻¹
                   fraction_of_calcite_not_dissolving_in_guts = (Z = 0.5, M = 0.75),
                   FeC_ratio_of_zooplankton = 10.0e-3,                                  #mmolFe molC⁻¹
                   FeZ_redfield_ratio = 3.0e-3,             #mmolFe molC⁻¹, remove this, is actually FeC_ratio_of_zooplankton
   
   
                   remineralisation_rate_of_DOC = 0.3 / day,                 #1/second
                   half_saturation_const_for_DOC_remin = 417.0,                #μmolCL⁻¹
                   NO3_half_saturation_const_for_DOC_remin = 0.03,            #μmolNL⁻¹
                   NH4_half_saturation_const_for_DOC_remin = 0.003,           #μmolNL⁻¹
                   PO4_half_saturation_const_for_DOC_remin = 0.003,       #μmolPL⁻¹
                   Fe_half_saturation_const_for_DOC_remin = 0.01,         #μmolFeL⁻¹
                   aggregation_rate_of_DOC_to_POC_1 = 0.37e-6 / day,          #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_DOC_to_POC_2 = 102.0e-6 / day,           #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_DOC_to_GOC_3 = 3530.0e-6 / day,          #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_DOC_to_POC_4 = 5095.0e-6 / day,          #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_DOC_to_POC_5 = 114.0e-6 / day,           #(μmolCL⁻¹)⁻¹s⁻¹
   
   
                   degradation_rate_of_POC = 0.025 / day,             #1/second
                   sinking_speed_of_POC = 2.0 / day,                    #ms⁻¹
                   min_sinking_speed_of_GOC = 30.0 / day,               #ms⁻¹
                   sinking_speed_of_dust = 2.0,                         #ms⁻¹
                   aggregation_rate_of_POC_to_GOC_6 = 25.9e-6 / day,     #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_POC_to_GOC_7 = 4452.0e-6 / day,     #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_POC_to_GOC_8 = 3.3e-6 / day,      #(μmolCL⁻¹)⁻¹s⁻¹
                   aggregation_rate_of_POC_to_GOC_9 = 47.1e-6 / day,     #(μmolCL⁻¹)⁻¹s⁻¹
                   min_scavenging_rate_of_iron = 3.0e-5 / day,          #1/second
                   slope_of_scavenging_rate_of_iron = 0.005 / day,    #d⁻¹μmol⁻¹L
                   scavenging_rate_of_iron_by_dust = 150.0 / day,       #s⁻¹mg⁻¹L
                   dissolution_rate_of_calcite = 0.197 / day,         #1/second
                   exponent_in_the_dissolution_rate_of_calcite = 1.0,
                   proportion_of_the_most_labile_phase_in_PSi = 0.5,
                   slow_dissolution_rate_of_PSi = 0.003 / day,        #1/second
                   fast_dissolution_rate_of_PSi = 0.025 / day,        #1/second
   
   
                   max_nitrification_rate = 0.05 / day,                           #1/sedonc
                   half_sat_const_for_denitrification1 = 1.0,                       #μmolO₂L⁻¹
                   half_sat_const_for_denitrification2 = 6.0,                       #μmolO₂L⁻¹
                   total_concentration_of_iron_ligands = 0.6,                     #nmolL⁻¹
                   max_rate_of_nitrogen_fixation = 0.013 / day,                         #μmolNL⁻¹s⁻¹
                   Fe_half_saturation_constant_of_nitrogen_fixation = 0.1,        #nmolFeL⁻¹
                   photosynthetic_parameter_of_nitrogen_fixation = 50.0,            #Wm⁻²
                   iron_concentration_in_sea_ice = 15.0,                            #nmolFeL⁻¹   
                   max_sediment_flux_of_Fe = 2.0 / day,                             #μmolFem⁻²s⁻¹
                   solubility_of_iron_in_dust = 0.02,
                   OC_for_ammonium_based_processes = 133/122,                     #molO₂(mol C)⁻¹
                   OC_ratio_of_nitrification = 32/122,                            #molO₂(mol C)⁻¹
                   CN_ratio_of_ammonification = 3/5,                              #molN(mol C)⁻¹
                   CN_ratio_of_denitrification = 105/16,                          #molN(mol C)⁻¹
                   NC_redfield_ratio = 16/122, 
                   PC_redfield_ratio = 1/122,                                   #molN(mol C)⁻¹
                   rain_ratio_parameter = 0.3,
                   bacterial_reference = 1.0,     #Not sure if this is what its called : denoted Bact_ref in paper

                   NC_stoichiometric_ratio_of_dentitrification = 0.86,
                   NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER = 0.0,     #again not sure what this is called
                   dissolution_rate_of_silicon = 1.0,
                   coefficient_of_bacterial_uptake_of_iron_in_POC = 0.5,
                   coefficient_of_bacterial_uptake_of_iron_in_GOC = 0.5,
                   max_FeC_ratio_of_bacteria = 10.0e-3,     #or 6
                   Fe_half_saturation_const_for_Bacteria = 0.03, #or 2.5e-10  

                   mixed_layer_shear = 1.0,
                   background_shear = 0.01, 

                   mixed_layer_depth = FunctionField{Center, Center, Center}(-100.0, grid),
                   euphotic_depth = Field{Center, Center, Nothing}(grid),
                   vertical_diffusivity  = ConstantField(1),
                   yearly_maximum_silicate = ConstantField(1),
                   dust_deposition = ZeroField(),

                  surface_photosynthetically_active_radiation = default_surface_PAR,

                  light_attenuation_model =
                    MultiBandPhotosyntheticallyActiveRadiation(; grid, 
                                                                 surface_PAR = surface_photosynthetically_active_radiation),

                  # just keep all this stuff for now but you can ignore it
                  sediment_model = nothing,

                  sinking_speeds = (POC = 2/day, 
                                    GOC = KernelFunctionOperation{Center, Center, Face}(DepthDependantSinkingSpeed(), 
                                                                                        grid, 
                                                                                        mixed_layer_depth, 
                                                                                        euphotic_depth)),
                  
                  carbonate_sat_ratio = ZeroField(),
                  open_bottom = true,

                  scale_negatives = false,

                  particles = nothing,
                  modifiers = nothing)

    if !isnothing(sediment_model) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    if (latitude isa Number) & !(grid isa RectilinearGrid)
        φ = φnodes(grid, Center(), Center(), Center())

        @warn "A latitude of $latitude was given but the grid has its own latitude ($(minimum(φ)), $(maximum(φ))) so the prescribed value is ignored"

        latitude = nothing 
    elseif isnothing(latitude) & (grid isa RectilinearGrid)
        throw(ArgumentError("You must prescribe a latitude when using a `RectilinearGrid`"))
    end

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
                                        Fe_half_saturation_const_for_Bacteria,    #not sure what this should be called

                                        mixed_layer_shear,
                                        background_shear,
 
                                        mixed_layer_depth,
                                        euphotic_depth,
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

@inline biogeochemical_auxiliary_fields(bgc::PISCES) = 
    (zₘₓₗ = bgc.mixed_layer_depth, 
     zₑᵤ = bgc.euphotic_depth, 
     Si̅ = bgc.yearly_maximum_silicate, 
     D_dust = bgc.dust_deposition, 
     Ω = bgc.carbonate_sat_ratio)

@inline required_biogeochemical_tracers(::PISCES) = 
    (:P, :D, :Z, :M, 
     :Pᶜʰˡ, :Dᶜʰˡ, 
     :Pᶠᵉ, :Dᶠᵉ, 
     :Dˢⁱ, 
     :DOC, :POC, :GOC, 
     :SFe, :BFe, 
     :PSi, 
     :NO₃, :NH₄, 
     :PO₄, :Fe, :Si, 
     :CaCO₃, :DIC, :Alk, 
     :O₂, 
     :T)

@inline required_biogeochemical_auxiliary_fields(::PISCES) =
    (:zₘₓₗ, :zₑᵤ, :Si̅, :D_dust, :Ω, :PAR, :PAR₁, :PAR₂, :PAR₃)

# for sinking things like POM this is how we tell oceananigans ther sinking speed
@inline function biogeochemical_drift_velocity(bgc::PISCES, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end

const small_particle_components = Union{Val{:POC}, Val{:SFe}}
const large_particle_components = Union{Val{:GOC}, Val{:BFe}, Val{:PSi}, Val{:CaCO₃}} 
# not sure what the point of PSi is since they don't use it to compute the sinking speed anymore anyway

biogeochemical_drift_velocity(bgc::PISCES, ::small_particle_components) = (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.POC)
biogeochemical_drift_velocity(bgc::PISCES, ::large_particle_components) = (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.GOC)

# don't worry about this for now
adapt_structure(to, pisces::PISCES) =
    PISCES(adapt(to, pisces.parameter_1),
           adapt(to, pisces.sinking_velocities))

# you can updatye these if you want it to have a pretty way of showing uyou its a pisces model
summary(::PISCES{FT}) where {FT} = string("PISCES{$FT}") 

show(io::IO, model::PISCES) = print(io, string("Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model")) # maybe add some more info here

@inline maximum_sinking_velocity(bgc::PISCES) = maximum(abs, bgc.sinking_velocities.bPOM.w) # might need ot update this for wghatever the fastest sinking pareticles are

include("phytoplankton.jl")
include("calcite.jl")
include("carbonate_system.jl")
include("dissolved_organic_matter.jl")
include("iron_in_particles.jl")
include("iron.jl")
include("nitrates_ammonium.jl")
include("oxygen.jl")
include("phosphates.jl")
include("particulate_organic_matter.jl")
include("silicon_in_particles.jl")
include("silicon.jl")
include("zooplankton.jl")

function update_biogeochemical_state!(model, bgc::PISCES)
    # this should come from utils
    #update_mixed_layer_depth!(bgc, model)

    PAR = biogeochemical_auxiliary_fields(model.biogeochemistry.light_attenuation).PAR

    compute_euphotic_depth!(bgc.euphotic_depth, PAR)

    # these should be model specific since they're not useful elsewhere
    #update_darkness_residence_time!(bgc, model)
    #update_yearly_maximum_silicate!(bgc, model)

    return nothing
end

# to work with the sediment model we need to tell in the redfield ratio etc. of some things, but for now we can ignore
@inline redfield(i, j, k, val_tracer_name, bgc::PISCES, tracers) = NaN

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline remineralisation_receiver(::PISCES) = :NH₄

# this is for positivity preservation, if you can work it out it would be great, I don't think PISCES conserves C but probably does Nitrogen
@inline conserved_tracers(::PISCES) = NaN

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here

@inline chlorophyll(bgc::PISCES, model) = model.tracers.Pᶜʰˡ + model.tracers.Dᶜʰˡ
end # module
