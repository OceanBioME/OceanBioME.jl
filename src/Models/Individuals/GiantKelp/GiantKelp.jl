module GiantKelpModel

export GiantKelp

using OceanBioME.Models.SugarKelpModel: LinearOptimalTemperatureRange

#=

     maximum_ammonia_uptake :: FT = 12 * 0.5 * 24 * 14 / (10^6)# 0.5 from structural_dry_weight_per_area
    ammonia_half_saturation :: FT = 1.3

     maximum_nitrate_uptake :: FT = 10 * 0.5 * 24 * 14 / (10^6)# 0.5 from structural_dry_weight_per_area
    nitrate_half_saturation :: FT = 4.0
    =#

@kwdef struct GiantKelp{GR, PP, RR, NU, ML, FT}
                          growth :: GR = ReserveAmmoniaGrowth(; maximum_specific_growth_rate = 0.2, # per macmods
                                                                low_area_enhancement_limit = Inf, 
                                                                base_growth_rate = 0.0, # get rid of area dependancy
                                                                seasonal_response = NotSeasonal(),
                                                                temperature_response = LinearOptimalTemperatureRange(; lower_optimal = 14.0,# per macmods
                                                                                                                       upper_optimal = 20.0,# per macmods
                                                                                                                       lower_gradient = 1/14,# per macmods
                                                                                                                       upper_gradient = -1/3))# per macmods
                  photosynthesis :: PP = ArrheniusPhotosynthesis() # can't constrain - I think the temperatures are likely to be wrong
                     respiration :: RR = BasalAndUptakeRespiration() # can't constrain
                 nitrogen_uptake :: NU = DroopUptake(; nitrate_half_saturation = 10.2, # per macmods
                                                       ammonia_half_saturation = 5.31, # per macmods - these seem high
                                                       maximum_nitrate_uptake = 752 / 10^6 * 14 / 100 / 0.0413, # per macmods -> gN / g / day
                                                       maximum_ammonia_uptake = 739 / 10^6 * 14 / 100 / 0.0413) # per macmods -> gN / g / day
                  mortality_loss :: ML = LinearMortality((0.0391+0.0199)/2) # mean from estimated parameters in SB canopy and water column https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.108
                                         # this is waaaaay larger than the sugar kelp parameterisation...SmallAreaLimitedMortality()
                                         # kind of matches macmods though - TODO: add wave dependance
        minimum_nitrogen_reserve :: FT = 0.01 # per macmods
        maximum_nitrogen_reserve :: FT = 0.04 # mer macmods
          minimum_carbon_reserve :: FT = 0.01 # arbitary guess

               structural_carbon :: FT = 0.19 # minimum in sbclter is 20%...
             structural_nitrogen :: FT = 0.0165 # sbclter minimum is 2.65% -> if min reserve is 1% per MACMODS we need this much in structure

  structural_dry_weight_per_area :: FT = 0.0413 #±0.006 per sblter data - quite the scatter though 
end

SugarKelpParticles(n; grid, kelp_parameters = NamedTuple(), kwargs...) = 
    BiogeochemicalParticles(n; grid, biogeochemistry = SugarKelp(; kelp_parameters...), kwargs...)

@inline required_particle_fields(::GiantKelp) = (:A, :N, :C)
@inline required_tracers(::GiantKelp) = (:u, :v, :w, :T, :NO₃, :NH₄, :PAR)
@inline coupled_tracers(::GiantKelp) = (:NO₃, :NH₄, :DIC, :O₂, :DOC, :DON, :bPOC, :bPON)

include("nitrogen_uptake.jl")
include("growth.jl")
include("seasonality.jl")
include("photosynthesis.jl")
include("respiration.jl")
include("mortality.jl")
include("equations.jl")

end # module