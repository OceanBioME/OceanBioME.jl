using Oceananigans.Grids: AbstractGrid
using Oceananigans.Units: day, days

# from https://web.mit.edu/globalchange/www/MITJPSPGC_Rpt122.pdf except nitrate and iron half saturations from 
# https://github.com/CliMA/ClimaOceanBiogeochemistry.jl/blob/main/src/carbon_alkalinity_nutrients.jl
# assumes all calcite goes into the particulate carbon, possibly won't conserve carbon when not used with variable redfield or the ParticulateDissolved below
@kwdef struct ImplicitProductivity{FT} # plankton
    maximum_community_productivity :: FT = 3/(365days) # mmol P / m³ / s
             light_half_saturation :: FT = 25.0        # W / m²
         phosphate_half_saturation :: FT = 0.5         # mmol P / m³
           nitrate_half_saturation :: FT = 7.17        # mmol N / m³
              iron_half_saturation :: FT = 1e-4        # mmol Fe / m³

       dissolved_fraction_of_waste :: FT = 0.67

                      carbon_ratio :: FT = 117.0       # mol C / mol P
                    nitrogen_ratio :: FT = 16.9        # mol N / mol P
                        iron_ratio :: FT = 4.68e-4     # mol Fe / mol P
                        rain_ratio :: FT = 0.07        # mol CaCO₃ / mol C
end

const ImplicitProductivityNPD = NutrientsPlanktonDetritus{<:Any, <:ImplicitProductivity}

required_biogeochemical_auxiliary_fields(::ImplicitProductivity) = (:PAR, )

@inline function community_productivity(bgc::ImplicitProductivityNPD, i, j, k, fields, auxiliary_fields)
    α  = bgc.plankton.maximum_community_productivity
    kₗ = bgc.plankton.light_half_saturation

    PAR = @inbounds auxiliary_fields.PAR[i, j, k]

    Lₙ = nutrient_limitation(bgc.nutrients, bgc.plankton, i, j, k, fields)
    Lₗ = PAR / (PAR + kₗ)

    return α * Lₙ * Lₗ
end

@inline phytoplankton_primary_production(bgc::ImplicitProductivityNPD, args...) =
    bgc.plankton.carbon_ratio * community_productivity(bgc, args...)

@inline solid_waste(bgc::ImplicitProductivityNPD, args...) =
    (1 - bgc.plankton.dissolved_fraction_of_waste) * community_productivity(bgc, args...)

@inline solid_carbon_waste(bgc::ImplicitProductivityNPD, args...) =
    (1 + bgc.plankton.rain_ratio) * bgc.plankton.carbon_ratio * solid_waste(bgc, args...)

@inline dissolved_waste(bgc::ImplicitProductivityNPD, args...) =
    (1 - bgc.plankton.dissolved_fraction_of_waste) * community_productivity(bgc, args...)

@inline plankton_organic_nitrogen_waste(bgc::ImplicitProductivityNPD, args...) = 
    bgc.plankton.nitrogen_ratio * dissolved_waste(bgc, args...)

@inline plankton_organic_carbon_waste(bgc::ImplicitProductivityNPD, args...) = 
    bgc.plankton.carbon_ratio * dissolved_waste(bgc, args...)

@inline calcite_uptake(bgc::ImplicitProductivityNPD, args...) =
    bgc.plankton.rain_ratio * (1 - bgc.plankton.dissolved_fraction_of_waste) * phytoplankton_primary_production(bgc, args...)

@inline nutrient_uptake(bgc::ImplicitProductivityNPD, i, j, k, val_name::Val{:N}, fields, auxiliary_fields) =
    community_productivity(bgc, i, j, k, fields, auxiliary_fields)

@inline nutrient_uptake(bgc::ImplicitProductivityNPD, i, j, k, val_name::Val{:NO₃}, fields, auxiliary_fields) =
    bgc.plankton.nitrogen_ratio * community_productivity(bgc, i, j, k, fields, auxiliary_fields)

@inline nutrient_uptake(bgc::ImplicitProductivityNPD, i, j, k, val_name::Val{:PO₄}, fields, auxiliary_fields) =
    community_productivity(bgc, i, j, k, fields, auxiliary_fields)  

@inline nutrient_uptake(bgc::ImplicitProductivityNPD, i, j, k, val_name::Val{:Fe}, fields, auxiliary_fields) =
    bgc.plankton.iron_ratio * community_productivity(bgc, i, j, k, fields, auxiliary_fields)

@inline grazing(bgc::ImplicitProductivityNPD, i, j, k, var_name, fields, auxiliary_fields) = 
    @inbounds zero(fields[1][1, 1, 1])

@inline plankton_inorganic_nitrogen_waste(bgc::ImplicitProductivityNPD, i, j, k, fields, auxiliary_fields) = 
    @inbounds zero(fields[1][1, 1, 1])

@inline plankton_inorganic_carbon_waste(bgc::ImplicitProductivityNPD, i, j, k, fields, auxiliary_fields) = 
    @inbounds zero(fields[1][1, 1, 1])

@inline plankton_inorganic_iron_waste(bgc::ImplicitProductivityNPD, i, j, k, fields, auxiliary_fields) = 
    @inbounds zero(fields[1][1, 1, 1])

@inline plankton_inorganic_phosphate_waste(bgc::ImplicitProductivityNPD, i, j, k, fields, auxiliary_fields) = 
    @inbounds zero(fields[1][1, 1, 1])

@inline function nutrient_limitation(::Nutrient, plankton, i, j, k, fields)
    k = plankton.phosphate_half_saturation
    N = @inbounds fields.N[i, j, k]

    return N/(N+k)
end

@inline function nutrient_limitation(::NitratePhosphateIron, plankton, i, j, k, fields)
    kₚ = plankton.phosphate_half_saturation
    kₙ = plankton.nitrate_half_saturation
    kᵢ = plankton.iron_half_saturation

    PO₄ = @inbounds fields.PO₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    Fe  = @inbounds  fields.Fe[i, j, k]

    return min(PO₄ / (PO₄ + kₚ), NO₃ / (NO₃ + kₙ), Fe / (Fe + kᵢ))
end

struct ParticulateDissolved{FT, SV} # detritus
      dissolved_remineralisation_rate :: FT
    particulate_remineralisation_rate :: FT
         particulate_sinking_velocity :: SV 

    function ParticulateDissolved(grid::AbstractGrid{FT};
                                  dissolved_remineralisation_rate::FT = 2/(365days),
                                  particulate_remineralisation_rate::FT = 5.7e-7,
                                  particulate_sinking_speed = 10/day,
                                  open_bottom = false) where FT
        particulate_sinking_velocity = setup_velocity_fields((; DOP = particulate_sinking_speed), grid, open_bottom; three_D = true).DOP

        SV = typeof(particulate_sinking_velocity)

        return new{FT, SV}(dissolved_remineralisation_rate,
                           particulate_remineralisation_rate,
                           particulate_sinking_velocity)
    end
end

required_biogeochemical_tracers(::ParticulateDissolved) = (:DOP, :POP)

biogeochemical_drift_velocity(pd::ParticulateDissolved, ::Val{:POP}) = 
    pd.particulate_sinking_velocity

const ParticulateDissolvedNPD =
    NutrientsPlanktonDetritus{<:Any, <:Any, <:ParticulateDissolved}

@inline (bgc::ParticulateDissolvedNPD)(i, j, k, grid, val_name::Val{:POP}, clock, fields, auxiliary_fields) = @inbounds (
    solid_waste(bgc, i, j, k, fields, auxiliary_fields) 
  - fields.POP[i, j, k] * bgc.detritus.particulate_remineralisation_rate
)

@inline (bgc::ParticulateDissolvedNPD)(i, j, k, grid, val_name::Val{:DOP}, clock, fields, auxiliary_fields) = @inbounds (
    dissolved_waste(bgc, i, j, k, fields, auxiliary_fields) 
  - fields.DOP[i, j, k] * bgc.detritus.dissolved_remineralisation_rate
)

@inline detritus_inorganic_waste(bgc::ParticulateDissolvedNPD, i, j, k, fields, auxiliary_fields) = @inbounds (
    fields.DOP[i, j, k] * bgc.detritus.dissolved_remineralisation_rate
  + fields.POP[i, j, k] * bgc.detritus.particulate_remineralisation_rate
)

const ImplicitProductivityParticulateDissolvedNPD =
    NutrientsPlanktonDetritus{<:Any, <:ImplicitProductivity, <:ParticulateDissolved}

@inline detritus_inorganic_carbon_waste(bgc::ImplicitProductivityParticulateDissolvedNPD, args...) =
    bgc.plankton.carbon_ratio * detritus_inorganic_waste(bgc, args...)

@inline detritus_inorganic_nitrogen_waste(bgc::ImplicitProductivityParticulateDissolvedNPD, args...) =
    bgc.plankton.carbon_ratio * detritus_inorganic_waste(bgc, args...)

@inline detritus_inorganic_phosphate_waste(bgc::ImplicitProductivityParticulateDissolvedNPD, args...) =
    detritus_inorganic_waste(bgc, args...)

@inline detritus_inorganic_iron_waste(bgc::ImplicitProductivityParticulateDissolvedNPD, args...) =
    bgc.plankton.iron_ratio * detritus_inorganic_waste(bgc, args...)

@inline calcite_dissolution(bgc::ImplicitProductivityParticulateDissolvedNPD, i, j, k, fields, auxiliary_fields) =
   @inbounds bgc.plankton.rain_ratio * bgc.plankton.carbon_ratio * fields.POP[i, j, k] * bgc.detritus.particulate_remineralisation_rate

# constructor

ImplicitBiology(grid;
                nutrients = NitratePhosphateIron(),
                plankton = ImplicitProductivity(),
                detritus = ParticulateDissolved(grid),
                carbonate_system = CarbonateSystem(),
                kwargs...) =
    NutrientsPlanktonDetritus(grid; 
                              nutrients,
                              plankton,
                              detritus,
                              carbonate_system,
                              kwargs...)

