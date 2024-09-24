using Oceananigans.Operators: ℑzᵃᵃᶜ

"""
    TwoCompartementCarbonIronParticles

A quota parameterisation for particulate organic matter with two size classes,
each with carbon and iron compartements, and a silicate compartement for the
large size class.

Confusingly we decided to name these compartmenets `POC` and `GOC` for the small
and large carbon classes, `SFe` and `BFe` for the small and  ̶l̶a̶r̶g̶e̶ big iron 
compartements, and `PSi` for the  ̶l̶a̶r̶g̶e̶ particulate silicon (*not* the 
phytoplankton silicon).
"""
@kwdef struct TwoCompartementCarbonIronParticles{FT, AP}
                      temperature_sensetivity :: FT = 1.066        #
                          base_breakdown_rate :: FT = 0.025 / day  # 1 / s
# (1 / (mmol C / m³),  1 / (mmol C / m³),  1 / (mmol C / m³) / s,  1 / (mmol C / m³) / s)
                       aggregation_parameters :: AP = (25.9, 4452, 3.3, 47.1) .* (10^-6 / day)
                 minimum_iron_scavenging_rate :: FT = 3e-5/day     # 1 / s
           load_specific_iron_scavenging_rate :: FT = 0.005/day    # 1 / (mmol C / m³) / s
  small_fraction_of_bacterially_consumed_iron :: FT = 0.5          #
  large_fraction_of_bacterially_consumed_iron :: FT = 0.5          #
                base_liable_silicate_fraction :: FT = 0.5          #
            fast_dissolution_rate_of_silicate :: FT = 0.025/day    # 1 / s
            slow_dissolution_rate_of_silicate :: FT = 0.003/day    # 1 / s

            # calcite
                              base_rain_ratio :: FT = 0.3         # 
                base_calcite_dissolution_rate :: FT = 0.197 / day # 1 / s
                 calcite_dissolution_exponent :: FT = 1.0         # 

            # iron in particles
               maximum_iron_ratio_in_bacteria :: FT = 0.06      # μmol Fe / mmol C
            iron_half_saturation_for_bacteria :: FT = 0.3       # μmol Fe / m³
                maximum_bacterial_growth_rate :: FT = 0.6 / day # 1 / s
end

const TwoCompartementPOCPISCES = PISCES{<:Any, <:Any, <:Any, TwoCompartementCarbonIronParticles}

required_biogeochemical_tracers(::TwoCompartementCarbonIronParticles) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃)

@inline edible_flux_rate(poc, i, j, k, grid, fields) = flux_rate(Val(:POC), i, j, k, grid, fields) + flux_rate(Val(:GOC), i, j, k, grid, fields)
@inline edible_iron_flux_rate(poc, i, j, k, grid, fields) = flux_rate(Val(:SFe), i, j, k, grid, fields) + flux_rate(Val(:BFe), i, j, k, grid, fields)

@inline flux_rate(::Val{:POC}, i, j, k, grid, fields) = @inbounds fields.POC[i, j, k] * ℑzᵃᵃᶜ(i, j, k, grid, fields.wPOC)
@inline flux_rate(::Val{:GOC}, i, j, k, grid, fields) = @inbounds fields.GOC[i, j, k] * ℑzᵃᵃᶜ(i, j, k, grid, fields.wGOC)
@inline flux_rate(::Val{:SFe}, i, j, k, grid, fields) = @inbounds fields.SFe[i, j, k] * ℑzᵃᵃᶜ(i, j, k, grid, fields.wPOC)
@inline flux_rate(::Val{:BFe}, i, j, k, grid, fields) = @inbounds fields.BFe[i, j, k] * ℑzᵃᵃᶜ(i, j, k, grid, fields.wGOC)

const small_particle_components = Union{Val{:POC}, Val{:SFe}}
const large_particle_components = Union{Val{:GOC}, Val{:BFe}, Val{:PSi}, Val{:CaCO₃}} 

biogeochemical_drift_velocity(bgc::TwoCompartementCarbonIronParticles, ::small_particle_components) = 
    (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.POC)

biogeochemical_drift_velocity(bgc::TwoCompartementCarbonIronParticles, ::large_particle_components) = 
    (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.GOC)

@inline function aggregation(poc::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields)
    a₁, a₂, a₃, a₄ = poc.aggregation_parameters

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear


    z = znode(i, j, k, grid, Center(), Center(), Center())

    zₘₓₗ = @inbounds fields.zₘₓₗ[i, j, k]
    
    POC = @inbounds fields.POC[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    return shear * (a₁ * POC^2 + a₂ * POC * GOC) + a₃ * POC * GOC + a₄ * POC^2
end

@inline function specific_degredation_rate(poc::TwoCompartementParticulateOrganicMatter, i, j, k, grid, bgc, clock, fields)
    λ₀ = poc.base_breakdown_rate
    b  = poc.temperature_sensetivity

    O₂ = @inbounds fields.O₂[i, j, k]

    ΔO₂ = anoxia_factor(bgc, O₂)

    return λ₀ * b^T * (1 - 0.45 * ΔO₂)
end

include("carbon.jl")
include("iron.jl")
include("silicate.jl")
include("calcite.jl")
include("nano_diatom_coupling.jl")
include("micro_meso_zoo_coupling.jl")