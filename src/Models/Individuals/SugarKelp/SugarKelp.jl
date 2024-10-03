module SugarKelpModel

export SugarKelp

using Roots

using Oceananigans.Units

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers

# disgsuting number of parameters
@kwdef struct SugarKelp{FT, TL}
    temperature_limit :: TL = LinearOptimalTemperatureRange() # TODO: split up more parameterisations like this
    growth_rate_adjustment :: FT = 4.5
    photosynthetic_efficiency :: FT = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)
    minimum_carbon_reserve :: FT = 0.01
    structural_carbon :: FT = 0.2
    exudation :: FT = 0.5
    erosion_exponent :: FT = 0.22
    base_erosion_rate :: FT = 10^-6
    saturation_irradiance :: FT = 90 * day/ (10 ^ 6)
    structural_dry_weight_per_area :: FT = 0.5
    structural_dry_to_wet_weight :: FT = 0.0785
    carbon_reserve_per_carbon :: FT = 2.1213
    nitrogen_reserve_per_nitrogen :: FT = 2.72
    minimum_nitrogen_reserve :: FT = 0.0126
    maximum_nitrogen_reserve :: FT = 0.0216
    growth_adjustment_2 :: FT = 0.039 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve))
    growth_adjustment_1 :: FT = 0.18 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)) - growth_adjustment_2
    maximum_specific_growth_rate :: FT = 0.18
    structural_nitrogen :: FT = 0.0146
    photosynthesis_at_ref_temp_1 :: FT = 1.22e-3 * 24
    photosynthesis_at_ref_temp_2 :: FT = 1.3e-3 * 24
    photosynthesis_ref_temp_1 :: FT = 285.0
    photosynthesis_ref_temp_2 :: FT = 288.0
    photoperiod_1 :: FT = 0.85
    photoperiod_2 :: FT = 0.3
    respiration_at_ref_temp_1 :: FT = 2.785e-4 * 24
    respiration_at_ref_temp_2 :: FT = 5.429e-4 * 24
    respiration_ref_temp_1 :: FT = 285.0
    respiration_ref_temp_2 :: FT = 290.0
    photosynthesis_arrhenius_temp :: FT = (1 / photosynthesis_ref_temp_1 - 1 / photosynthesis_ref_temp_2) ^ -1 * log(photosynthesis_at_ref_temp_2 / photosynthesis_at_ref_temp_1)
    photosynthesis_low_temp :: FT = 271.0
    photosynthesis_high_temp :: FT = 296.0
    photosynthesis_high_arrhenius_temp :: FT = 1414.87
    photosynthesis_low_arrhenius_temp :: FT = 4547.89
    respiration_arrhenius_temp :: FT = (1 / respiration_ref_temp_1 - 1 / respiration_ref_temp_2) ^ -1 * log(respiration_at_ref_temp_2 / respiration_at_ref_temp_1)
    current_speed_for_0p65_uptake :: FT = 0.03
    nitrate_half_saturation :: FT = 4.0
    ammonia_half_saturation :: FT = 1.3
    maximum_nitrate_uptake :: FT = 10 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    maximum_ammonia_uptake :: FT = 12 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    current_1 :: FT = 0.72
    current_2 :: FT = 0.28
    current_3 :: FT = 0.045
    base_activity_respiration_rate :: FT = 1.11e-4 * 24
    base_basal_respiration_rate :: FT = 5.57e-5 * 24
    exudation_redfield_ratio :: FT = Inf
    adapted_latitude :: FT = 57.5
end 

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