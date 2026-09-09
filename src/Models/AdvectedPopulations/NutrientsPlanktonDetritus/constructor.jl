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
"""
    MARBL(grid; kwargs...)

The MARBL biogeochemistry, in its CESM2.1 configuration: three phytoplankton functional types (small
phytoplankton, diatoms and diazotrophs) and one zooplankton, with variable cell quotas in carbon,
chlorophyll, phosphorus, iron and silicon, a full grazing matrix, dissolved organic matter with
refractory pools, the dissolved iron and ligand system, and oxygen with nitrification and
denitrification.

Particulates sink the way MARBL sinks them: through the implicit mineral-protection (ballast) flux of
Armstrong et al. (2000). Particulate organic carbon and phosphorus, biogenic silica, particulate iron,
dust and calcite are therefore **not tracers** — each step the column is swept from the surface down
and the resulting remineralisation is what the tendencies read. MARBL has no explicit sinking speed.

This is a preset assembling the general components — [`ManyPhytoZoo`](@ref),
[`MultiElementRefractoryDissolved`](@ref) with a [`Ballast`](@ref), [`ImplicitExplicitCalcite`](@ref)
(the calcite-tracer-free sibling of [`ExplicitCalciumCarbonate`](@ref), since the sweep solves calcite),
[`ComplexedIron`](@ref) and [`RedoxOxygen`](@ref) — with the CESM2.1 parameter values. Any of them can be
replaced through the corresponding keyword.

!!! warning "Sibling components"
    [`MultiElementRefractoryDissolved`](@ref) pairs with [`ImplicitExplicitCalcite`](@ref), and
    [`MultiElementRefractoryDissolvedParticulate`](@ref) with [`ExplicitCalciumCarbonate`](@ref).
    If you replace one of the pair, replace the other: mixing them is a silent physics error rather
    than a load error.

Temperature is *not* a biogeochemical tracer: supply it from the physics, or as `T` in the model's
tracers.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `open_bottom`: whether the ballast floor flux leaves the domain (for a sediment to catch); forwarded
  to the [`Ballast`](@ref)
- `ballast`: the [`Ballast`](@ref) sinking parameters, which also carry the surface dust flux and the
  sedimentary iron flux
- any other keyword argument is forwarded to [`NutrientsPlanktonDetritus`](@ref)

See [`MARBL_Cocco`](@ref) for the four-PFT variant with coccolithophores,
[`MARBL_ExplicitSinking`](@ref) for the explicitly sinking testing configuration, and
[`BurialDenitrificationSediment`](@ref) to close the bottom boundary.
"""
MARBL(grid::AbstractGrid{FT};
      open_bottom = true,
      ballast = Ballast(FT; open_bottom),
      nutrients = Nutrients(NitrateAmmonia{FT}(; nitrification_rate = zero(FT)),  # nitrification via `RedoxOxygen`
                            PO₄, ComplexedIron(FT), Si),
      plankton = ManyPhytoZoo(FT),
      detritus = MultiElementRefractoryDissolved(grid; ballast),
      inorganic_carbon = ImplicitExplicitCalcite(grid),
      oxygen = RedoxOxygen(FT),
      surface_PAR = default_surface_PAR,
      light_attenuation = default_light(grid, surface_PAR),
      kwargs...) where FT =
    NutrientsPlanktonDetritus(grid; nutrients, plankton, detritus, inorganic_carbon, oxygen,
                              light_attenuation, kwargs...)

"""
    MARBL_ExplicitSinking(grid; sinking_speed = 10.0, kwargs...)

The MARBL biogeochemistry with the **explicitly sinking** particulate detritus
([`MultiElementRefractoryDissolvedParticulate`](@ref)) and an explicit calcite tracer
([`ExplicitCalciumCarbonate`](@ref)) in place of the implicit ballast flux of [`MARBL`](@ref).
Everything else — the plankton, nutrients, iron and ligand system, and oxygen chemistry — is as
[`MARBL`](@ref).

!!! warning "This is a testing configuration, not CESM2.1"
    MARBL has no explicit sinking: in MARBL sinking *is* the implicit ballast flux, so use
    [`MARBL`](@ref) for the CESM2.1 configuration. This preset exists because carrying the
    particulates as tracers makes a closed box conserve every element exactly, which is what the
    conservation tests need. Its `sinking_speed` and the particulate remineralisation rates of
    [`MultiElementRefractoryDissolvedParticulate`](@ref) are **provisional placeholders, not CESM2.1
    values**, and it has never been compared against the MARBL Fortran driver.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `sinking_speed`: the (provisional) sinking speed of every particulate class (m/day)
- `open_bottom`: whether particulates sink out of the bottom of the domain
- any other keyword argument is forwarded to [`MARBL`](@ref)
"""
MARBL_ExplicitSinking(grid::AbstractGrid{FT};
                      open_bottom = true,
                      sinking_speed = 10.0,   # m/day, a placeholder — MARBL has no explicit sinking
                      detritus = MultiElementRefractoryDissolvedParticulate(grid; sinking_speed, open_bottom),
                      inorganic_carbon = ExplicitCalciumCarbonate(grid),
                      kwargs...) where FT =
    MARBL(grid; open_bottom, detritus, inorganic_carbon, kwargs...)

"""
    MARBL_Cocco(grid; kwargs...)

The MARBL biogeochemistry in its `+cocco` configuration: as [`MARBL`](@ref), but with a fourth
phytoplankton functional type, coccolithophores, which are the sole calcifier and the only carbon-
limited type. The small phytoplankton do not calcify in this configuration, and several shared
parameters take different values.

Coccolithophore growth and calcification depend on aqueous CO₂, so the plankton allocates a field for it
and recomputes it each step from the carbonate system.
"""
MARBL_Cocco(grid::AbstractGrid{FT};
            plankton = MARBL_cocco_plankton(grid, FT),
            kwargs...) where FT =
    MARBL(grid; plankton, kwargs...)
