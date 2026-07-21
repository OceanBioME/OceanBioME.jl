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

for (nutrient, symbol) in pairs((phosphate = :PO₄, iron = :Fe, silicate = :Si))
    fname = Symbol(nutrient, :_limitation)
    @eval begin
        @inline function $fname(i, j, k, grid,
                                ::SingleTracerNutrient, 
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
