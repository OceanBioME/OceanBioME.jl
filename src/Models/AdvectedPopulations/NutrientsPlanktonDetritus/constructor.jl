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
                              negative_tracers = nothing,
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
- `negative_tracers`: treatment for negative tracer values, such as [`IgnoreNegativeTracerValues`](@ref),
  [`ClipNegativeTracers`](@ref), or [`ScaleNegativeTracers`](@ref)
- `scale_negatives`: convenience for `negative_tracers = ScaleNegativeTracers(; invalid_fill_value)`
- `invalid_fill_value`: value used by the `scale_negatives` convenience
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
                                   negative_tracers = nothing,
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
        isnothing(negative_tracers) ||
            throw(ArgumentError("Specify either `scale_negatives=true` or `negative_tracers`, not both."))
        negative_tracers = ScaleNegativeTracers(; invalid_fill_value)
    end

    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation,
                           sediment,
                           particles,
                           modifiers,
                           negative_tracers)
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

"""
    MITgcmDIC(grid; open_bottom = true, kwargs...)

Construct a [`NutrientsPlanktonDetritus`](@ref) preset matching the MITgcm DIC package
(Dutkiewicz et al., 2005). It tracks phosphate (`PO₄`), dissolved organic phosphorus (`DOP`),
particulate organic phosphorus (`POP`), iron (`Fe`) with ligand equilibrium and scavenging
([`SimpleIron`](@ref)), dissolved inorganic carbon (`DIC`), alkalinity (`Alk`), and oxygen (`O₂`).

Community productivity is computed by [`ImplicitProductivity`](@ref) (equivalent to MITgcm's
`bio_export.F`), limited by phosphate, iron, and light.

External iron sources (aeolian dust deposition, sediment flux) should be applied as Oceananigans
`Forcing` on the `Fe` tracer.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `open_bottom`: whether detritus can sink out of the bottom of the domain
- `surface_PAR`: the surface photosynthetically active radiation passed to the default light model
- any other keyword argument is forwarded to [`NutrientsPlanktonDetritus`](@ref) (e.g. `nutrients`,
  `plankton`, `detritus`, `oxygen`, `light_attenuation`)
"""
MITgcmDIC(grid::AbstractGrid{FT};
          open_bottom = true,
          nutrients = Nutrients(nothing, PO₄, SimpleIron{FT}(), nothing),
          plankton = ImplicitProductivity(FT;
                                          maximum_community_productivity = 2 / (360 * day),  # mmol P / m³ / s
                                          light_half_saturation = 30.0,                      # W / m²
                                          dissolved_fraction_of_waste = 0.67,
                                          carbon_ratio = 117.0,                              # mol C / mol P
                                          nitrogen_ratio = 16.0,                             # mol N / mol P
                                          iron_ratio = 4.68e-4,                              # mol Fe / mol P
                                          rain_ratio = 0.07,                                 # mol CaCO₃ / mol C
                                          nutrient_half_saturations =
                                              (phosphate = 0.5,                              # mmol P / m³
                                               iron = 1.2e-4)),                              # mmol Fe / m³
          detritus = DissolvedParticulate(grid, :DOP, :POP;
                                          dissolved_remineralisation_rate = 1 / (6 * 30 * day), # 1/s
                                          particulate_remineralisation_rate = 0.03 / day,
                                          dissolved_fraction_of_remineralisation = 0.0,
                                          sinking_speeds = 10 / day,
                                          open_bottom),
          inorganic_carbon = CarbonateSystem(),
          oxygen = Oxygen(FT;
                          production_oxygen_carbon_ratio = 170 / 117,     # |R_OP / R_CP|
                          nitrification_oxygen_carbon_ratio = 16 / 117),  # R_NP / R_CP
          surface_PAR = default_surface_PAR,
          light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR),
          kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, inorganic_carbon, oxygen, light_attenuation, kwargs...)