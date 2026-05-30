using Oceananigans.Units

using .PlanktonModels: limiting_nutrients

function NutrientsPlanktonDetritus(grid::AbstractGrid{FT}; 
                                   nutrients = Nutrients(nothing, nothing, nothing, nothing),
                                   plankton = Abiotic(),
                                   detritus = InstantRemineralisation(),
                                   inorganic_carbon = nothing,
                                   oxygen = nothing,
                                   light_attenuation = nothing,
                                   sediment = nothing,
                                   scale_negatives = false,
                                   invalid_fill_value = convert(FT, NaN),
                                   particles = nothing,
                                   modifiers = nothing) where FT

    underlying_biogeochemistry = 
        NutrientsPlanktonDetritus{eltype(grid)}(nutrients, 
                                                plankton, 
                                                detritus, 
                                                inorganic_carbon, 
                                                oxygen)

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
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end

ImplicitBiology(grid::AbstractGrid{FT};
                nutrients = Nutrients(N, PO₄, Fe, nothing),
                plankton = ImplicitProductivity(FT),
                detritus = DissolvedParticulate(grid, :DOP, :POP;
                                                dissolved_remineralisation_rate = 2/365/day,
                                                particulate_remineralisation_rate = 0.03/day,
                                                dissolved_fraction_of_remineralisation = 0.0,
                                                sinking_speeds = 10/day),
                inorganic_carbon = CarbonateSystem(),
                kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, inorganic_carbon, kwargs...)

NPZD(grid::AbstractGrid{FT};
     limiting_nutrients = (:nitrate, :ammonia),
     nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia(FT) : N, 
                          :phosphate in limiting_nutrients ? PO₄ : nothing, 
                          :iron in limiting_nutrients ? Fe : nothing, 
                          nothing),
     plankton = PhytoZoo(FT;
                        nutrient_half_saturations = (nitrate = 0.7,                     # mmol N/m³
                                                     ammonia = 0.001,                   # mmol N/m³
                                                     iron = 2e-4)[limiting_nutrients]), # mmol Fe / m³
     detritus = Detritus(grid),
     kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, kwargs...)

LOBSTER(grid::AbstractGrid{FT};
        limiting_nutrients = (:nitrate, :ammonia),
        nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia(FT) : N, 
                              :phosphate in limiting_nutrients ? PO₄ : nothing, 
                              :iron in limiting_nutrients ? Fe : nothing, 
                              nothing),
        plankton = PhytoZoo(FT;
                            nutrient_half_saturations = (nitrate = 0.7,                     # mmol N/m³
                                                         ammonia = 0.001,                   # mmol N/m³
                                                         iron = 2e-4)[limiting_nutrients]), # mmol Fe / m³
        detritus = DissolvedParticulate(grid),
        kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, kwargs...)