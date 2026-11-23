"""
    ImplicitProductivity

An implicit plankton which computes a community productivity
limited by available nutrients and light, but not tracking any
planktonic biomass.

Mainly from MITGCM https://web.mit.edu/globalchange/www/MITJPSPGC_Rpt122.pdf

Units: mmol P / m³
Tracers: none
Default limiting nutrients: N or NO₃, P, Fe

The `iron_ratio` (mol Fe / mol P) is either a number, the Fe:P of all the organic matter, or an [`IronQuota`](@ref).
The `rain_ratio` (mol CaCO₃ / mol C) is either a number or a [`GrowthDependentRainRatio`](@ref).
"""
struct ImplicitProductivity{FT, LN, TT, IR, RR} <: AbstractPlankton{LN}
    maximum_community_productivity :: FT
             light_half_saturation :: FT
       dissolved_fraction_of_waste :: FT
                      carbon_ratio :: FT
                    nitrogen_ratio :: FT
                        iron_ratio :: IR
                        rain_ratio :: RR

         nutrient_half_saturations :: NamedTuple{LN, TT}
end

"""
    IronQuota(FT = Float64; maximum_ratio = 35e-6 * 117, minimum_ratio = 2.7e-6 * 117, saturating_iron = 3e-4)

An Fe:P (mol Fe / mol P) for [`ImplicitProductivity`](@ref) which varies with the dissolved iron `Fe` (mmol Fe / m³),
as MARBL's Fe:C for growth (`gQfe`):

    Q = max(minimum_ratio, min(maximum_ratio, maximum_ratio × Fe / saturating_iron))

The defaults are MARBL's small phytoplankton's `gQfe_0` and `gQfe_min` at a C:P of 117, and `gQ_Fe_kFe_thres` × `kFe`.

The iron of the solid waste, Q times it, goes into `SimpleIron`'s particulate iron (`PFe`). The dissolved and
particulate organic matter carry no iron, so the iron of the dissolved waste is never taken up.
"""
struct IronQuota{FT}
      maximum_ratio :: FT # mol Fe / mol P
      minimum_ratio :: FT # mol Fe / mol P
    saturating_iron :: FT # mmol Fe / m³
end

IronQuota(FT = Float64; maximum_ratio = 35e-6 * 117, minimum_ratio = 2.7e-6 * 117, saturating_iron = 3e-4) =
    IronQuota(convert(FT, maximum_ratio), convert(FT, minimum_ratio), convert(FT, saturating_iron))

convert_iron_ratio(FT, ratio::Number) = convert(FT, ratio)
convert_iron_ratio(FT, q::IronQuota) =
    IronQuota(FT; maximum_ratio = q.maximum_ratio, minimum_ratio = q.minimum_ratio, saturating_iron = q.saturating_iron)

"""
    GrowthDependentRainRatio(FT = Float64;
                             maximum_ratio = 0.095,
                             growth_grazing_ratio = 15,
                             minimum_temperature = -2,
                             saturating_temperature = 4)

A calcium carbonate rain ratio (mol CaCO₃ / mol C) for [`ImplicitProductivity`](@ref) which varies with the
community's growth limitation `L` and the temperature `T` (°C):

    R = maximum_ratio × clamp((T - minimum_temperature) / (saturating_temperature - minimum_temperature), 0, 1) / (1 + (growth_grazing_ratio × L)²)

`L` is the community productivity over its `maximum_community_productivity`, so the light × nutrient limitation.
`1 / (1 + (growth_grazing_ratio × L)²)` is the small phytoplankton's share of the production after Dunne et al. (2005),
as in BLING (Galbraith et al. 2010), where the large phytoplankton's is (μ / λ)² / (1 + (μ / λ)²) for growth and
grazing rates with μ / λ = `growth_grazing_ratio` × `L`. The temperature term is that of MARBL's implicit calcifier,
with its `CaCO3_temp_thres2` and `CaCO3_temp_thres1`.

It needs the temperature tracer `T`, and calcium carbonate which sinks on its own (e.g. `ImplicitCalciumCarbonate`).
"""
struct GrowthDependentRainRatio{FT}
             maximum_ratio :: FT # mol CaCO₃ / mol C
      growth_grazing_ratio :: FT # –
       minimum_temperature :: FT # °C
    saturating_temperature :: FT # °C
end

GrowthDependentRainRatio(FT = Float64; maximum_ratio = 0.095, growth_grazing_ratio = 15,
                         minimum_temperature = -2, saturating_temperature = 4) =
    GrowthDependentRainRatio(convert(FT, maximum_ratio), convert(FT, growth_grazing_ratio),
                             convert(FT, minimum_temperature), convert(FT, saturating_temperature))

convert_rain_ratio(FT, ratio::Number) = convert(FT, ratio)
convert_rain_ratio(FT, r::GrowthDependentRainRatio) =
    GrowthDependentRainRatio(FT; maximum_ratio = r.maximum_ratio, growth_grazing_ratio = r.growth_grazing_ratio,
                             minimum_temperature = r.minimum_temperature, saturating_temperature = r.saturating_temperature)

function ImplicitProductivity(FT = Float64; 
                              maximum_community_productivity = 3/(365days), # mmol P / m³ / s
                              light_half_saturation = 25.0,                 # W / m²
                              dissolved_fraction_of_waste = 0.67,
                              carbon_ratio = 106.0,                         # mol C / mol P
                              nitrogen_ratio = 16.0,                        # mol N / mol P
                              iron_ratio = 4.68e-4,                         # mol Fe / mol P
                              rain_ratio = 0.07,                            # mol CaCO₃ / mol C

                              nutrient_half_saturations = 
                                  (phosphate = 0.5,                        # mmol P / m³
                                   nitrate = 7.17,                         # mmol N / m³
                                   iron = 1e-4))                           # mmol Fe / m³)

    limiting_nutrients = keys(nutrient_half_saturations)

    nutrient_half_saturations = NamedTuple{limiting_nutrients}(map(val->convert(FT, val), values(nutrient_half_saturations)))

    return ImplicitProductivity(convert(FT, maximum_community_productivity),
                                convert(FT, light_half_saturation),
                                convert(FT, dissolved_fraction_of_waste),
                                convert(FT, carbon_ratio),
                                convert(FT, nitrogen_ratio),
                                convert_iron_ratio(FT, iron_ratio),
                                convert_rain_ratio(FT, rain_ratio),
                                nutrient_half_saturations)
end

required_biogeochemical_tracers(::ImplicitProductivity) = tuple()
required_biogeochemical_auxiliary_fields(::ImplicitProductivity) = (:PAR, )

@inline phosphate_ratio(::ImplicitProductivity, ::NPD{FT}) where FT = one(FT)
@inline nitrogen_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.nitrogen_ratio
@inline carbon_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.carbon_ratio
@inline iron_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.iron_ratio
@inline iron_ratio(::ImplicitProductivity{<:Any, <:Any, <:Any, <:IronQuota}, ::NPD{FT}) where FT = zero(FT)
@inline calcium_carbonate_rain_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.rain_ratio

const GrowthDependentRainRatioProductivity = ImplicitProductivity{<:Any, <:Any, <:Any, <:Any, <:GrowthDependentRainRatio}

required_biogeochemical_tracers(::GrowthDependentRainRatioProductivity) = (:T, )

@inline function growth_dependent_rain_ratio(r::GrowthDependentRainRatio, L, T)
    temperature_limitation = clamp((T - r.minimum_temperature) / (r.saturating_temperature - r.minimum_temperature), 0, 1)
    small_plankton_share = 1 / (1 + (r.growth_grazing_ratio * L)^2)

    return r.maximum_ratio * temperature_limitation * small_plankton_share
end

@inline function calcium_carbonate_rain_ratio(i, j, k, grid, plankton::GrowthDependentRainRatioProductivity, bgc, fields, auxiliary_fields)
    T = @inbounds fields.T[i, j, k]
    L = community_productivity(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) / plankton.maximum_community_productivity

    return growth_dependent_rain_ratio(plankton.rain_ratio, L, T)
end

@inline iron_quota(q::IronQuota, Fe) = max(q.minimum_ratio, min(q.maximum_ratio, q.maximum_ratio * Fe / q.saturating_iron))

const IronQuotaProductivity = ImplicitProductivity{<:Any, <:Any, <:Any, <:IronQuota}

@inline gross_iron_uptake(i, j, k, grid, plankton::IronQuotaProductivity, bgc::NPD, fields, auxiliary_fields) =
    iron_quota(plankton.iron_ratio, @inbounds(fields.Fe[i, j, k])) *
    community_productivity(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

@inline function particulate_iron_production(i, j, k, grid, plankton::IronQuotaProductivity, bgc::NPD, fields, auxiliary_fields)
    Q = iron_quota(plankton.iron_ratio, @inbounds(fields.Fe[i, j, k]))

    return Q * solid_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
end

@inline limiting_nutrients(::ImplicitProductivity{<:Any, <:NamedTuple{LN}}) where LN = LN

@inline nutrient_half_saturations(implicit::ImplicitProductivity, ::Val{:N})   = implicit.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(implicit::ImplicitProductivity, ::Val{:NO₃}) = implicit.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(implicit::ImplicitProductivity, ::Val{:NH₄}) = implicit.nutrient_half_saturations.ammonia
@inline nutrient_half_saturations(implicit::ImplicitProductivity, ::Val{:PO₄}) = implicit.nutrient_half_saturations.phosphate
@inline nutrient_half_saturations(implicit::ImplicitProductivity, ::Val{:Fe})  = implicit.nutrient_half_saturations.iron

@inline inorganic_waste(i, j, k, grid, ::ImplicitProductivity, bgc::NPD{FT}, args...) where FT = zero(FT)

@inline dissolved_waste(i, j, k, grid, plankton::ImplicitProductivity, bgc, args...) =
    plankton.dissolved_fraction_of_waste * community_productivity(i, j, k, grid, plankton, bgc, args...)

@inline solid_waste(i, j, k, grid, plankton::ImplicitProductivity, bgc::NPD{FT}, args...) where FT =
    (one(FT) - plankton.dissolved_fraction_of_waste) * community_productivity(i, j, k, grid, plankton, bgc, args...)

@inline nutrient_uptake(i, j, k, grid, plankton::ImplicitProductivity, bgc, fields, auxiliary_fields) =
    community_productivity(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) 

@inline function community_productivity(i, j, k, grid, plankton::ImplicitProductivity, bgc, fields, auxiliary_fields)
    α  = plankton.maximum_community_productivity
    kₗ = plankton.light_half_saturation

    PAR = @inbounds auxiliary_fields.PAR[i, j, k]

    Lₙ = nutrient_limitation(i, j, k, grid, bgc.nutrients, plankton, bgc, fields, auxiliary_fields)
    Lₗ = PAR / (PAR + kₗ)

    return α * Lₙ * Lₗ
end

@inline nutrient_uptake(i, j, k, grid, 
                        ::Val{:NO₃},
                        plankton::ImplicitProductivity, bgc, 
                        fields, auxiliary_fields) =
    nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) * nitrogen_ratio(i, j, k, grid, plankton, bgc, fields)

@inline nutrient_uptake(i, j, k, grid, 
                        ::Val{:NH₄},
                        ::ImplicitProductivity, ::NPD{FT}, 
                        fields, auxiliary_fields) where FT = zero(FT)

@inline calcium_carbonate_precipitation(i, j, k, grid,
                              bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:ImplicitProductivity},
                              fields, auxiliary_fields) = (
    solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)

@inline chlorophyll(::ImplicitProductivity, model) = ZeroField()

# admin

Adapt.adapt_structure(to, ip::ImplicitProductivity) = 
    ImplicitProductivity(
        adapt(to, ip.maximum_community_productivity),
        adapt(to, ip.light_half_saturation),
        adapt(to, ip.dissolved_fraction_of_waste),
        adapt(to, ip.carbon_ratio),
        adapt(to, ip.nitrogen_ratio),
        adapt(to, ip.iron_ratio),
        adapt(to, ip.rain_ratio),
        adapt(to, ip.nutrient_half_saturations)
    )

Base.summary(::ImplicitProductivity) = string("ImplicitProductivity")

function Base.show(io::IO, ip::ImplicitProductivity)
    msg = summary(ip) * "\n"
    msg *= "├── Community productivity (≤ $(round(ip.maximum_community_productivity*365*day, digits = 1)) mmol P/year) straight to detritus\n"
    msg *= "├── C:N:P:Fe = $(ip.carbon_ratio):$(ip.nitrogen_ratio):1:$(ip.iron_ratio)\n"
    msg *= "└── Limiting nutrients: $(limiting_nutrients(ip))"

    print(io, msg)

    return nothing
end
