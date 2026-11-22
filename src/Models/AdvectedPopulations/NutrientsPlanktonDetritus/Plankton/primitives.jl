# building blocks

abstract type AbstractPlankton{LN} end

# growth limitation - to use these primitives plankton be AbstractPlankton with LN = limiting nutrients
# a tuple of limiting nutrients
# and `nutrient_half_saturations(plankton, symbol)`
@inline nutrient_limitation(i, j, k, grid, nutrients, args...) = 
    min(nitrogen_limitation(i, j, k, grid, nutrients.nitrogen, args...),
        phosphate_limitation(i, j, k, grid, nutrients.phosphate, args...),
        iron_limitation(i, j, k, grid, nutrients.iron, args...),
        silicate_limitation(i, j, k, grid, nutrients.silicate, args...))

const IronNutrient = Union{SingleTracerNutrient, SimpleIron}

for (nutrient, symbol, dispatch) in ((:phosphate, :PO₄, :SingleTracerNutrient),
                                     (     :iron,  :Fe,  :IronNutrient),
                                     ( :silicate,  :Si,  :SingleTracerNutrient))
    fname = Symbol(nutrient, :_limitation)
    @eval begin
        @inline function $fname(i, j, k, grid,
                                ::$dispatch,
                                plankton::AbstractPlankton{LN},
                                bgc::NutrientsPlanktonDetritus{FT},
                                fields, auxiliary_fields) where {LN, FT}
            if $(QuoteNode(nutrient)) in LN
                kN = nutrient_half_saturations(plankton, Val($(QuoteNode(symbol))))
                N = fields.$symbol[i, j, k]

                return N / (N + kN)
            else
                one(FT)
            end
        end

        $fname(i, j, k, grid,
               ::Nothing,
               plankton,
               ::NutrientsPlanktonDetritus{FT},
               fields, auxiliary_fields) where FT =
            one(FT)
    end
end

@inline function nitrogen_limitation(i, j, k, grid,
                                     ::SingleTracerNutrient,
                                     plankton::AbstractPlankton{LN},
                                     ::NutrientsPlanktonDetritus{FT},
                                     fields, auxiliary_fields) where {LN, FT}

    if (:nitrate in LN) | (:nitrogen in LN)
        kN = nutrient_half_saturations(plankton, Val(:N)) 
        N = fields.N[i, j, k]

        return N / (N + kN)
    else
        one(FT)
    end
end

@inline nitrogen_limitation(i, j, k, grid,
                            ::Nothing, 
                            plankton,
                            ::NutrientsPlanktonDetritus{FT},
                            fields, auxiliary_fields) where FT = one(FT)

@inline function nitrogen_limitation(i, j, k, grid,
                                     ::NitrateAmmonia, 
                                     plankton::AbstractPlankton{LN},
                                     ::NutrientsPlanktonDetritus{FT},
                                     fields, auxiliary_fields) where {LN, FT}

    if (:nitrate in LN) & (:ammonia in LN)
        kNO₃ = nutrient_half_saturations(plankton, Val(:NO₃)) 
        kNH₄ = nutrient_half_saturations(plankton, Val(:NH₄))

        @inbounds begin 
            NO₃ = fields.NO₃[i, j, k]
            NH₄ = fields.NH₄[i, j, k]
        end

        return (NO₃/kNO₃ + NH₄/kNH₄) / (1 + NO₃/kNO₃ + NH₄/kNH₄)
    elseif :nitrogen in LN
        kN = nutrient_half_saturations(plankton, Val(:NO₃)) 
        N = fields.NO₃[i, j, k]

        return N / (N + kN)
    else
        return one(FT)
    end
end

struct MondoLightLimitation end
struct AnalyticalLightLimitation end # This formulation is justified because you can rearrange to αI/√(μ² + α²I²) which if you put simple I(t) in can be integrated analytically

@inline light_limitation(::MondoLightLimitation, PAR, kPAR) = PAR / (kPAR + PAR)
@inline light_limitation(::AnalyticalLightLimitation, PAR, kPAR) = PAR / sqrt(PAR^2 + kPAR^2)

# temperature limitation - `temperature_dependence` is `nothing` (no dependence), a bare Q₁₀ number
# (legacy `PhytoZoo` behaviour), or any object callable on the temperature in °C

@inline temperature_limitation(::Nothing, i, j, k, fields) = 1

@inline function temperature_limitation(Q10::Number, i, j, k, fields)
    T = @inbounds fields.T[i, j, k]

    return Q10 ^ (T / 10)
end

@inline function temperature_limitation(response, i, j, k, fields)
    T = @inbounds fields.T[i, j, k]

    return response(T)
end

"""
    EppleyCurve(FT = Float64; temperature_coefficient = 1.88, reference_temperature = 0)

An exponential (Eppley, 1972) temperature response, `Q₁₀^((T - Tᵣ) / 10)`, where `Q₁₀` is the
`temperature_coefficient` and `Tᵣ` the `reference_temperature` (°C). The default `Q₁₀` of 1.88 is
Eppley's 1.066ᵀ curve, which approximately doubles growth every 10 °C.
"""
struct EppleyCurve{FT}
    temperature_coefficient :: FT
      reference_temperature :: FT
end

EppleyCurve(FT = Float64; temperature_coefficient = 1.88, reference_temperature = 0) =
    EppleyCurve(convert(FT, temperature_coefficient), convert(FT, reference_temperature))

@inline (curve::EppleyCurve)(T) = curve.temperature_coefficient ^ ((T - curve.reference_temperature) / 10)

convert_temperature_dependence(FT, ::Nothing) = nothing
convert_temperature_dependence(FT, curve::EppleyCurve) = EppleyCurve(FT; curve.temperature_coefficient, curve.reference_temperature)
convert_temperature_dependence(FT, response) = response

temperature_dependence_summary(::Nothing) = "none"
temperature_dependence_summary(response) = summary(response)

Base.summary(curve::EppleyCurve) = "EppleyCurve(Q₁₀ = $(curve.temperature_coefficient), Tᵣ = $(curve.reference_temperature)°C)"
Base.show(io::IO, curve::EppleyCurve) = print(io, summary(curve))
