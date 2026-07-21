using Oceananigans.Units

using OceanBioME.Light: 
    TwoBandPhotosyntheticallyActiveRadiation,
    PrescribedAttenuationPAR

using Oceananigans.Fields: ConstantField

using .PlanktonModels: limiting_nutrients

"""
    NutrientsPlanktonDetritus(grid;
                              nutrients = Nutrients(nothing, nothing, nothing, nothing),
                              plankton = Abiotic(),
                              detritus = InstantRemineralisationDetritus(),
                              inorganic_carbon = nothing,
                              oxygen = nothing,
                              light_attenuation = nothing,
                              sediment = nothing,
                              scale_negatives = false,
                              invalid_fill_value = NaN,
                              particles = nothing,
                              modifiers = nothing)

Construct a biogeochemical model in the modular Nutrients-Plankton-Detritus (NPD) framework by
assembling it from pluggable components. Each slot may be swapped independently, and the preset
[`LOBSTER`](@ref), [`NPZD`](@ref), and [`ImplicitBiology`](@ref)
constructors are just particular choices of these components.

The set of tracers the model evolves is determined by the components you choose.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to configure sinking speeds
- `nutrients`: the inorganic nutrient pool(s) that limit growth; a [`Nutrients`](@ref) grouping of
  nitrogen, phosphate, iron, and silicate components. Defaults to no explicit nutrients (all slots
  `nothing`), in which case nutrients are implicitly conserved and not tracked
- `plankton`: the planktonic (living) component, defaults to [`Abiotic`](@ref) (no biology). Options
  include [`PhytoZoo`](@ref) and [`ImplicitProductivity`](@ref)
- `detritus`: the non-living organic component, defaults to [`InstantRemineralisationDetritus`](@ref)
  which returns waste straight to the nutrient pool. Options include [`Detritus`](@ref),
  [`DissolvedParticulate`](@ref), and [`CarbonNitrogenDissolvedParticulate`](@ref)
- `inorganic_carbon`: optional inorganic carbon component, defaults to `nothing`; can be a
  [`CarbonateSystem`](@ref) which adds dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`)
- `oxygen`: optional oxygen component, defaults to `nothing`; can be an [`Oxygen`](@ref) which adds
  oxygen (`O₂`)
- `light_attenuation`: light attenuation model which integrates the attenuation of available light
- `sediment`: slot for a sediment model (`AbstractSediment`)
- `scale_negatives`: whether to add a [`ScaleNegativeTracers`](@ref) modifier to keep tracers non-negative
- `invalid_fill_value`: the value used to fill invalid tracer values when `scale_negatives` is `true`
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry after the tendencies have been
  calculated or when the state is updated
"""
function NutrientsPlanktonDetritus(grid::AbstractGrid{FT};
                                   nutrients = Nutrients(nothing, nothing, nothing, nothing),
                                   plankton = Abiotic(),
                                   detritus = InstantRemineralisationDetritus(),
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

"""
    ImplicitBiology(grid; limiting_nutrients = (:nitrate, :iron, :phosphate), open_bottom = true, kwargs...)

Construct the `ImplicitBiology` preset of the [`NutrientsPlanktonDetritus`](@ref) framework: a model
which computes community productivity limited by available nutrients and light without explicitly
tracking planktonic biomass (see [`ImplicitProductivity`](@ref)). It is suited to large-scale or
long-timescale simulations where resolving full plankton dynamics is too expensive.

By default it couples the [`Nutrients`](@ref) selected by `limiting_nutrients`, an
[`ImplicitProductivity`](@ref) plankton, a two-class [`DissolvedParticulate`](@ref) detritus
(`DOP`/`POP`), a [`CarbonateSystem`](@ref), and a [`PrescribedAttenuationPAR`](@ref) light model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `limiting_nutrients`: a tuple of the nutrients that limit productivity, drawn from `:nitrate` (or
  `:ammonia` to split nitrogen into `NO₃`/`NH₄`), `:phosphate`, and `:iron`
- `open_bottom`: whether detritus can sink out of the bottom of the domain
- `surface_PAR`: the surface photosynthetically active radiation passed to the default light model
- any other keyword argument is forwarded to [`NutrientsPlanktonDetritus`](@ref) (e.g. `oxygen`,
  `nutrients`, `plankton`, `detritus`, `light_attenuation`)
"""
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
                                                     phosphate = 0.5,                    # mmol N/m³
                                                     iron = 1e-4)[limiting_nutrients]),  # mmol Fe / m³),
                detritus = DissolvedParticulate(grid, :DOP, :POP;
                                                dissolved_remineralisation_rate = 2/365/day,
                                                particulate_remineralisation_rate = 0.03/day,
                                                dissolved_fraction_of_remineralisation = 0.0,
                                                sinking_speeds = 10/day,
                                                open_bottom),
                inorganic_carbon = CarbonateSystem(),
                surface_PAR = default_surface_PAR,
                light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR),
                kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, inorganic_carbon, light_attenuation, kwargs...)

"""
    NPZD(grid; limiting_nutrients = (:nitrate,), open_bottom = true, kwargs...)

Construct the `NPZD` (Nutrient-Phytoplankton-Zooplankton-Detritus) preset of the
[`NutrientsPlanktonDetritus`](@ref) framework. It couples the [`Nutrients`](@ref) selected by
`limiting_nutrients` with a [`PhytoZoo`](@ref) plankton (phytoplankton `P` and zooplankton `Z`)
parameterised after Kuhn et al. (2015), a single-class [`Detritus`](@ref) pool (`D`), and the default
two-band light model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `limiting_nutrients`: a tuple of the nutrients that limit growth, drawn from `:nitrate` (or
  `:ammonia` to split nitrogen into `NO₃`/`NH₄`), `:phosphate`, and `:iron`
- `open_bottom`: whether detritus can sink out of the bottom of the domain
- `surface_PAR`: the surface photosynthetically active radiation passed to the default light model
- any other keyword argument is forwarded to [`NutrientsPlanktonDetritus`](@ref)
"""
NPZD(grid::AbstractGrid{FT};
     limiting_nutrients = (:nitrate, ),
     open_bottom = true,
     nutrients = Nutrients(:ammonia in limiting_nutrients ? NitrateAmmonia{FT}() : N, 
                           :phosphate in limiting_nutrients ? PO₄ : nothing, 
                           :iron in limiting_nutrients ? Fe : nothing, 
                           nothing),
     plankton = PhytoZoo(grid;
                         nutrient_half_saturations = (nitrate = 2.3868,                     # mmol N/m³
                                                      ammonia = 0.001,                   # mmol N/m³
                                                      iron = 2e-4)[limiting_nutrients], # mmol Fe / m³
                         phytoplankton_maximum_growth_rate = 0.6989 / day,
                         phytoplankton_exudation_fraction = zero(FT),
                         temperature_coefficient = 1.88,
                         phytoplankton_mortality_rate = (0.066 + 0.0101)/day,
                         preference_for_phytoplankton = one(FT),
                         grazing_half_saturation = 0.5573,
                         zooplankton_mortality_rate = 0.3395 / day,
                         zooplankton_excretion_rate = 0.0102 / day,
                         zooplankton_assimilation_fraction = 0.9116,
                         phytoplankton_sinking_speed = 0.2551/day,
                         excretion_inorganic_fraction = one(FT),
                         phytoplankton_solid_waste_fraction = 0.0101 / (0.066 + 0.0101),
                         maximum_grazing_rate = 2.1522 / day,
                         light_limitation = PlanktonModels.AnalyticalLightLimitation(),
                         light_half_saturation = (0.6989/day)/(0.1953/day)),
     detritus = Detritus(grid; open_bottom),
     surface_PAR = default_surface_PAR,
     light_attenuation = default_light(grid, surface_PAR),
     kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, light_attenuation, kwargs...)

"""
    LOBSTER(grid; limiting_nutrients = (:nitrate, :ammonia), open_bottom = true, kwargs...)

Construct the [LOBSTER](@ref) preset of the [`NutrientsPlanktonDetritus`](@ref) framework. By
default it splits nitrogen into nitrate (`NO₃`) and ammonia (`NH₄`) via [`NitrateAmmonia`](@ref),
couples a [`PhytoZoo`](@ref) plankton (phytoplankton `P` and zooplankton `Z`), a dissolved-and-two-
particulate [`DissolvedParticulate`](@ref) detritus (`DOM`, `sPOM`, `bPOM`), and the default two-band
light model. Pass `inorganic_carbon = CarbonateSystem()` and/or `oxygen = Oxygen()` to add carbonate
chemistry and oxygen.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `limiting_nutrients`: a tuple of the nutrients that limit growth, drawn from `:nitrate`, `:ammonia`,
  `:phosphate`, and `:iron`
- `open_bottom`: whether detritus can sink out of the bottom of the domain
- `surface_PAR`: the surface photosynthetically active radiation passed to the default light model
- any other keyword argument is forwarded to [`NutrientsPlanktonDetritus`](@ref) (e.g.
  `inorganic_carbon`, `oxygen`)
"""
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
        light_attenuation = default_light(grid, surface_PAR),
        kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, light_attenuation, kwargs...)