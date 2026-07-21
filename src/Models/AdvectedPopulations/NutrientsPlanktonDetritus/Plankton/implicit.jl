"""
    ImplicitProductivity

An implicit plankton which computes a community productivity
limited by available nutrients and light, but not tracking any
planktonic biomass.

Mainly from MITGCM https://web.mit.edu/globalchange/www/MITJPSPGC_Rpt122.pdf

Units: mmol P / m³
Tracers: none
Default limiting nutrients: N or NO₃, P, Fe
"""
struct ImplicitProductivity{FT, LN, TT} <: AbstractPlankton{LN} 
    maximum_community_productivity :: FT
             light_half_saturation :: FT
       dissolved_fraction_of_waste :: FT
                      carbon_ratio :: FT
                    nitrogen_ratio :: FT
                        iron_ratio :: FT
                        rain_ratio :: FT

         nutrient_half_saturations :: NamedTuple{LN, TT}
end

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
                                convert(FT, iron_ratio),
                                convert(FT, rain_ratio),
                                nutrient_half_saturations)
end

required_biogeochemical_tracers(::ImplicitProductivity) = tuple()
required_biogeochemical_auxiliary_fields(::ImplicitProductivity) = (:PAR, )

@inline phosphate_ratio(::ImplicitProductivity, ::NPD{FT}) where FT = one(FT)
@inline nitrogen_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.nitrogen_ratio
@inline carbon_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.carbon_ratio
@inline iron_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.iron_ratio
@inline calcite_rain_ratio(plankton::ImplicitProductivity, ::NPD{FT}) where FT = plankton.rain_ratio

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

@inline calcite_precipitation(i, j, k, grid,
                              bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:ImplicitProductivity},
                              fields, auxiliary_fields) = (
    solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)

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
