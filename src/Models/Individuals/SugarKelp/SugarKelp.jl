"""
Sugar kelp model of [Broch2012](@citet) and updated by [Broch2013](@citet), [Fossberg2018](@citet), and [Broch2019](@citet).

Prognostic properties
=====================
* Area: A (dm²)
* Nitrogen reserve: N (gN/gSW)
* Carbon reserve: C (gC/gSW)

Tracer dependencies
===================
* Nitrates: NO₃ (mmol N/m³)
* Ammonia: NH₄ (mmol N/m³)
* Photosynthetically available radiation: PAR (einstein/m²/day)
* Temperature: T (°C)

""" 
module SugarKelpModel

export SugarKelp, SugarKelpParticles

using Oceananigans.Grids: AbstractGrid
using Oceananigans.Units

using OceanBioME: NewtonRaphsonSolver
using OceanBioME.Particles: BiogeochemicalParticles

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers

"""
    SugarKelp

Defines the parameters for `SugarKelp` biogeochemistry.
"""
struct SugarKelp{FT, TL, SO}
                     temperature_limit :: TL
                growth_rate_adjustment :: FT
             photosynthetic_efficiency :: FT
                minimum_carbon_reserve :: FT
                     structural_carbon :: FT
                             exudation :: FT
                      erosion_exponent :: FT
                     base_erosion_rate :: FT
                 saturation_irradiance :: FT
        structural_dry_weight_per_area :: FT
          structural_dry_to_wet_weight :: FT
             carbon_reserve_per_carbon :: FT
         nitrogen_reserve_per_nitrogen :: FT
              minimum_nitrogen_reserve :: FT
              maximum_nitrogen_reserve :: FT
                   growth_adjustment_2 :: FT
                   growth_adjustment_1 :: FT
          maximum_specific_growth_rate :: FT
                   structural_nitrogen :: FT
          photosynthesis_at_ref_temp_1 :: FT
          photosynthesis_at_ref_temp_2 :: FT 
             photosynthesis_ref_temp_1 :: FT
             photosynthesis_ref_temp_2 :: FT
                         photoperiod_1 :: FT
                         photoperiod_2 :: FT
             respiration_at_ref_temp_1 :: FT
             respiration_at_ref_temp_2 :: FT
                respiration_ref_temp_1 :: FT
                respiration_ref_temp_2 :: FT
         photosynthesis_arrhenius_temp :: FT
               photosynthesis_low_temp :: FT
              photosynthesis_high_temp :: FT
    photosynthesis_high_arrhenius_temp :: FT
     photosynthesis_low_arrhenius_temp :: FT
            respiration_arrhenius_temp :: FT
         current_speed_for_0p65_uptake :: FT
               nitrate_half_saturation :: FT
               ammonia_half_saturation :: FT
                maximum_nitrate_uptake :: FT
                maximum_ammonia_uptake :: FT
                             current_1 :: FT
                             current_2 :: FT
                             current_3 :: FT 
        base_activity_respiration_rate :: FT
           base_basal_respiration_rate :: FT
              exudation_redfield_ratio :: FT
                      adapted_latitude :: FT
                                solver :: SO

    function SugarKelp(FT = Float64; 
                       temperature_limit::TL = LinearOptimalTemperatureRange{FT}(),
                       growth_rate_adjustment = 4.5,
                       photosynthetic_efficiency = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60),
                       minimum_carbon_reserve = 0.01,
                       structural_carbon = 0.2,
                       exudation = 0.5,
                       erosion_exponent = 0.22,
                       base_erosion_rate = 10^-6,
                       saturation_irradiance = 90 * day/ (10 ^ 6),
                       structural_dry_weight_per_area = 0.5,
                       structural_dry_to_wet_weight = 0.0785,
                       carbon_reserve_per_carbon = 2.1213,
                       nitrogen_reserve_per_nitrogen = 2.72,
                       minimum_nitrogen_reserve = 0.0126,
                       maximum_nitrogen_reserve = 0.0216,
                       growth_adjustment_2 = 0.039 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)),
                       growth_adjustment_1 = 0.18 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)) - growth_adjustment_2,
                       maximum_specific_growth_rate = 0.18,
                       structural_nitrogen = 0.0146,
                       photosynthesis_at_ref_temp_1 = 1.22e-3 * 24,
                       photosynthesis_at_ref_temp_2 = 1.3e-3 * 24,
                       photosynthesis_ref_temp_1 = 285.0,
                       photosynthesis_ref_temp_2 = 288.0,
                       photoperiod_1 = 0.85,
                       photoperiod_2 = 0.3,
                       respiration_at_ref_temp_1 = 2.785e-4 * 24,
                       respiration_at_ref_temp_2 = 5.429e-4 * 24,
                       respiration_ref_temp_1 = 285.0,
                       respiration_ref_temp_2 = 290.0,
                       photosynthesis_arrhenius_temp = (1 / photosynthesis_ref_temp_1 - 1 / photosynthesis_ref_temp_2) ^ -1 * log(photosynthesis_at_ref_temp_2 / photosynthesis_at_ref_temp_1),
                       photosynthesis_low_temp = 271.0,
                       photosynthesis_high_temp = 296.0,
                       photosynthesis_high_arrhenius_temp = 1414.87,
                       photosynthesis_low_arrhenius_temp = 4547.89,
                       respiration_arrhenius_temp = (1 / respiration_ref_temp_1 - 1 / respiration_ref_temp_2) ^ -1 * log(respiration_at_ref_temp_2 / respiration_at_ref_temp_1),
                       current_speed_for_0p65_uptake = 0.03,
                       nitrate_half_saturation = 4.0,
                       ammonia_half_saturation = 1.3,
                       maximum_nitrate_uptake = 10 / structural_dry_weight_per_area * 24 * 14 / (10^6),
                       maximum_ammonia_uptake = 12 / structural_dry_weight_per_area * 24 * 14 / (10^6),
                       current_1 = 0.72,
                       current_2 = 0.28,
                       current_3 = 0.045,
                       base_activity_respiration_rate = 1.11e-4 * 24,
                       base_basal_respiration_rate = 5.57e-5 * 24,
                       exudation_redfield_ratio = Inf,
                       adapted_latitude = 57.5,
                       solver::SO = NewtonRaphsonSolver{FT, Int}(; atol = eps(FT(1e-9)))) where {TL, SO}
    
        @warn "The sugar kelp model includes an inate seasonality which is formulated assuming that the model time starts at the start of the year"

        return new{FT, TL, SO}(temperature_limit, 
                               growth_rate_adjustment, photosynthetic_efficiency,
                               minimum_carbon_reserve, structural_carbon,
                               exudation, erosion_exponent, base_erosion_rate,
                               saturation_irradiance,
                               structural_dry_weight_per_area, structural_dry_to_wet_weight,
                               carbon_reserve_per_carbon, nitrogen_reserve_per_nitrogen,
                               minimum_nitrogen_reserve, maximum_nitrogen_reserve,
                               growth_adjustment_2, growth_adjustment_1, maximum_specific_growth_rate,
                               structural_nitrogen,
                               photosynthesis_at_ref_temp_1, photosynthesis_at_ref_temp_2,
                               photosynthesis_ref_temp_1, photosynthesis_ref_temp_2, photoperiod_1, photoperiod_2,
                               respiration_at_ref_temp_1, respiration_at_ref_temp_2, respiration_ref_temp_1, respiration_ref_temp_2,
                               photosynthesis_arrhenius_temp, photosynthesis_low_temp, photosynthesis_high_temp,
                               photosynthesis_high_arrhenius_temp, photosynthesis_low_arrhenius_temp,
                               respiration_arrhenius_temp,
                               current_speed_for_0p65_uptake, nitrate_half_saturation, ammonia_half_saturation,
                               maximum_nitrate_uptake, maximum_ammonia_uptake, current_1, current_2, current_3, 
                               base_activity_respiration_rate, base_basal_respiration_rate, exudation_redfield_ratio,
                               adapted_latitude, 
                               solver)
    end
end 

"""
    SugarKelpParticles(n; grid, kelp_parameters = NamedTuple(), kwargs...)

Sets up `n` sugar kelp `BiogeochemicalParticles` with default parameters except
those specified in `kelp_parameters`. `kwagrs` are passed onto `BiogeochemicalParticles`.
"""
SugarKelpParticles(n; grid::AbstractGrid{FT}, kelp_parameters = NamedTuple(), kwargs...) where FT = 
    BiogeochemicalParticles(n; grid, biogeochemistry = SugarKelp(FT; kelp_parameters...), kwargs...)

@inline required_particle_fields(::SugarKelp) = (:A, :N, :C)
@inline required_tracers(::SugarKelp) = (:u, :v, :w, :T, :NO₃, :NH₄, :PAR)
@inline coupled_tracers(::SugarKelp) = (:NO₃, :NH₄, :DIC, :O₂, :DOC, :DON, :bPOC, :bPON)
# can overload like:
# @inline coupled_tracers(::SugarKelp) = (NO₃ = :NO₃, NH₄ = :NH₄, DON = :DOM, bPON = :POM)
# if you have a simpler model ... there must be a better way todo this

include("equations.jl")
include("coupling.jl")
include("show.jl")

end # module
