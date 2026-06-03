using Oceananigans.Units

using OceanBioME.Light: 
    TwoBandPhotosyntheticallyActiveRadiation,
    PrescribedPhotosyntheticallyActiveRadiation

using Oceananigans.Fields: ConstantField

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
        scaler = ScaleNegativeTracers(underlying_biogeochemistry; invalid_fill_value)
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

const default_light = TwoBandPhotosyntheticallyActiveRadiation
const default_surface_PAR = 100

ImplicitBiology(grid::AbstractGrid{FT};
                limiting_nutrients = (:nitrate, :iron, :phosphate),
                open_bottom = true,
                nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia{FT}() : N, 
                                      :phosphate in limiting_nutrients ? PO₄ : nothing, 
                                      :iron in limiting_nutrients ? Fe : nothing, 
                                      nothing),
                plankton = ImplicitProductivity(FT;
                                                nutrient_half_saturations = 
                                                    (nitrate = 7.17,                     # mmol N/m³
                                                     phosphate = 0.5,                   # mmol N/m³
                                                     iron = 1e-4)[limiting_nutrients]), # mmol Fe / m³),
                detritus = DissolvedParticulate(grid, :DOP, :POP;
                                                dissolved_remineralisation_rate = 2/365/day,
                                                particulate_remineralisation_rate = 0.03/day,
                                                dissolved_fraction_of_remineralisation = 0.0,
                                                sinking_speeds = 10/day,
                                                open_bottom),
                inorganic_carbon = CarbonateSystem(),
                surface_PAR = default_surface_PAR,
                light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(surface_PAR)),
                kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, inorganic_carbon, light_attenuation, kwargs...)

NPZD(grid::AbstractGrid{FT};
     limiting_nutrients = (:nitrate, ),
     open_bottom = true,
     nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia{FT}() : N, 
                           :phosphate in limiting_nutrients ? PO₄ : nothing, 
                           :iron in limiting_nutrients ? Fe : nothing, 
                           nothing),
     plankton = PhytoZoo(FT;
                         nutrient_half_saturations = (nitrate = 0.7,                     # mmol N/m³
                                                      ammonia = 0.001,                   # mmol N/m³
                                                      iron = 2e-4)[limiting_nutrients]), # mmol Fe / m³
     detritus = Detritus(grid; open_bottom),
     surface_PAR = default_surface_PAR,
     light_attenuation = default_light(; grid, surface_PAR),
     kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, light_attenuation, kwargs...)

LOBSTER(grid::AbstractGrid{FT};
        limiting_nutrients = (:nitrate, :ammonia),
        open_bottom = true,
        nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia{FT}() : N, 
                              :phosphate in limiting_nutrients ? PO₄ : nothing, 
                              :iron in limiting_nutrients ? Fe : nothing, 
                              nothing),
        plankton = PhytoZoo(FT;
                            nutrient_half_saturations = (nitrate = 0.7,                     # mmol N/m³
                                                         ammonia = 0.001,                   # mmol N/m³
                                                         iron = 2e-4)[limiting_nutrients]), # mmol Fe / m³
        detritus = DissolvedParticulate(grid; open_bottom),
        surface_PAR = default_surface_PAR,
        light_attenuation = default_light(; grid, surface_PAR),
        kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, light_attenuation, kwargs...)