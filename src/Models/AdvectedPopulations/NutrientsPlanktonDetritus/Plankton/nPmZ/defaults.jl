"""
    Autotroph([FT = Float64;] kwargs...)

One phytoplankton functional type for a [`ManyPhytoZoo`](@ref) plankton component. Its process
*traits* — whether it silicifies, calcifies (implicitly or explicitly), fixes nitrogen, is carbon
limited, and which temperature form it uses — are type parameters, so every branch that depends on
them is resolved at compile time and the PFT's tracer set follows from them: `<name>C`, `<name>Chl`,
`<name>P`, `<name>Fe`, plus `<name>Si` for a silicifier and `<name>CaCO₃` for a calcifier.

The same name may carry different traits in different configurations, so traits live per-instance
rather than in a name table.

Keyword Arguments
=================

- `silicifier`, `implicit_calcifier`, `explicit_calcifier`, `carbon_limited`, `nitrogen_fixer`: process traits
- `temperature_function`: `:q10` or `:power` (Eppley)
- `nitrate`, `ammonia`, `phosphate`, `DOP`, `iron`, `silicate`: nutrient half-saturations (mmol/m³)
- `photosynthesis_rate_reference_per_day`, `initial_PI_slope_per_day`,
  `maximum_chlorophyll_to_nitrogen`: growth and photoacclimation parameters
- `temperature_sensitivity`, `temperature_reference`, `temperature_threshold`: temperature response
- `linear_mortality_rate_per_day`, `quadratic_mortality_rate_per_day`,
  `maximum_aggregation_rate_per_day`, `minimum_aggregation_rate_per_day`,
  `particulate_loss_fraction`, `loss_threshold_concentration`,
  `cold_loss_threshold_concentration`: loss parameters
- `growth_iron_quota_reference`, `minimum_iron_quota`: iron quota parameters
- `carbon_dioxide_half_saturation`: only used when `carbon_limited`
"""
function Autotroph(FT = Float64;
                   silicifier           = false,
                   implicit_calcifier   = false,
                   explicit_calcifier   = false,
                   carbon_limited       = false,
                   nitrogen_fixer       = false,
                   temperature_function = :q10,

                   nitrate, ammonia, phosphate, DOP, iron, silicate,   # mmol/m³

                   photosynthesis_rate_reference_per_day,              # 1/day
                   initial_PI_slope_per_day = 0.39,                    # mmol C m²/(mg Chl W day)
                   maximum_chlorophyll_to_nitrogen,                    # mg Chl/mmol N

                   temperature_sensitivity = 1.7,                      # Q10
                   temperature_reference   = 30.0,                     # °C
                   temperature_threshold,                              # °C

                   linear_mortality_rate_per_day    = 0.1,             # 1/day
                   quadratic_mortality_rate_per_day = 0.01,            # 1/day / (mmol C/m³)
                   maximum_aggregation_rate_per_day = 0.5,             # 1/day
                   minimum_aggregation_rate_per_day,                   # 1/day
                   particulate_loss_fraction = 0.0,

                   loss_threshold_concentration,                       # mmol C/m³
                   cold_loss_threshold_concentration,                  # mmol C/m³

                   growth_iron_quota_reference,                        # mmol Fe/mmol C
                   minimum_iron_quota,                                 # mmol Fe/mmol C

                   carbon_dioxide_half_saturation = 0.0)               # mmol/m³, 0 unless carbon limited

    traits = (; silicifier, implicit_calcifier, explicit_calcifier, carbon_limited, nitrogen_fixer, temperature_function)

    khs = (nitrate   = convert(FT, nitrate),
           ammonia   = convert(FT, ammonia),
           phosphate = convert(FT, phosphate),
           DOP       = convert(FT, DOP),
           iron      = convert(FT, iron),
           silicate  = convert(FT, silicate))

    return Autotroph{traits, typeof(khs), FT}(
        khs,
        convert(FT, photosynthesis_rate_reference_per_day / day),
        convert(FT, initial_PI_slope_per_day / day),
        convert(FT, maximum_chlorophyll_to_nitrogen),
        convert(FT, temperature_sensitivity),
        convert(FT, temperature_reference),
        convert(FT, temperature_threshold),
        convert(FT, linear_mortality_rate_per_day / day),
        convert(FT, quadratic_mortality_rate_per_day / day),
        convert(FT, maximum_aggregation_rate_per_day / day),
        convert(FT, minimum_aggregation_rate_per_day / day),
        convert(FT, particulate_loss_fraction),
        convert(FT, loss_threshold_concentration),
        convert(FT, cold_loss_threshold_concentration),
        convert(FT, growth_iron_quota_reference),
        convert(FT, minimum_iron_quota),
        convert(FT, carbon_dioxide_half_saturation))
end

"""
    Grazing([FT = Float64;] prey, kwargs...)

One cell of the grazing matrix: a predator's relationship with ONE prey class. `prey` lists the class
members (a `Symbol` for the usual single-member class, or a tuple), which may be autotrophs and/or
zooplankton. The class is grazed as a single pool and the flux is split back over its members in
proportion to their biomass.

Keyword Arguments
=================

- `prey`: the member name(s) of this prey class
- `maximum_grazing_rate_per_day`, `grazing_half_saturation`: the grazing functional response
- `sigmoidal`: `true` for a sigmoidal response, `false` (default) for Michaelis–Menten
- `fraction_to_zooplankton`, `fraction_to_particulate`, `fraction_to_dissolved`: routing of the
  grazed carbon (the remainder goes to dissolved inorganic carbon)
- `detritus_fraction`: the food-weighted share of this predator's own losses routed to particulates
"""
function Grazing(FT = Float64;
                 prey,
                 maximum_grazing_rate_per_day,      # 1/day
                 grazing_half_saturation = 1.2,     # mmol C/m³
                 fraction_to_zooplankton,
                 fraction_to_particulate,
                 fraction_to_dissolved = 0.06,
                 detritus_fraction,
                 sigmoidal = false)

    prey = prey isa Symbol ? (prey, ) : Tuple(prey)

    return Grazing{prey, FT}(convert(FT, maximum_grazing_rate_per_day / day),
                             convert(FT, grazing_half_saturation),
                             convert(FT, fraction_to_zooplankton),
                             convert(FT, fraction_to_particulate),
                             convert(FT, fraction_to_dissolved),
                             convert(FT, detritus_fraction),
                             sigmoidal)
end

"""
    Heterotroph([FT = Float64;] grazing, kwargs...)

One zooplankton for a [`ManyPhytoZoo`](@ref) plankton component, carrying carbon biomass (tracer
`<name>C`) and its row of the grazing matrix.

The phosphorus and iron quotas are *not* here: they are shared by every zooplankton and live on
[`ManyPhytoZoo`](@ref), because once a zooplankter eats another a per-predator quota would leave an
unrouted P/Fe surplus and break conservation.

Keyword Arguments
=================

- `grazing`: a `NamedTuple` of [`Grazing`](@ref), one per prey class this predator grazes
- `linear_mortality_rate_per_day`, `quadratic_mortality_rate_per_day`,
  `loss_threshold_concentration`: loss parameters
- `temperature_sensitivity`, `temperature_reference`: the predator's own Q10 response
"""
function Heterotroph(FT = Float64;
                     grazing,
                     linear_mortality_rate_per_day    = 0.1,    # 1/day
                     quadratic_mortality_rate_per_day = 0.4,    # 1/day / (mmol C/m³)
                     loss_threshold_concentration     = 0.075,  # mmol C/m³
                     temperature_sensitivity          = 1.7,    # Q10
                     temperature_reference            = 30.0)   # °C

    return Heterotroph{typeof(grazing), FT}(
        grazing,
        convert(FT, linear_mortality_rate_per_day / day),
        convert(FT, quadratic_mortality_rate_per_day / day),
        convert(FT, loss_threshold_concentration),
        convert(FT, temperature_sensitivity),
        convert(FT, temperature_reference))
end

"""
    ManyPhytoZoo([FT = Float64;] kwargs...)

A multi-PFT plankton component for the `plankton` slot of a [`NutrientsPlanktonDetritus`](@ref)
model, carrying an arbitrary set of phytoplankton ([`Autotroph`](@ref)) and zooplankton
([`Heterotroph`](@ref)) with a full grazing matrix between them.

The currency is **carbon**: each autotroph carries carbon, chlorophyll, phosphorus and iron (plus
silicon and calcite where its traits call for them) as separate tracers with variable cell quotas,
and each zooplankton carries carbon alone. The tracer names are built from the PFT names, so any
naming works.

Parameters shared by every PFT live here; per-PFT parameters live on the [`Autotroph`](@ref) and
[`Heterotroph`](@ref) themselves. Defaults are the CESM2.1 values; the [`MARBL`](@ref) preset
assembles the standard PFT sets.

Keyword Arguments
=================

- `autotrophs`, `zooplankton`: `NamedTuple`s of [`Autotroph`](@ref)/[`Heterotroph`](@ref), keyed by
  the PFT names that give the tracers their prefixes
- `carbon_dioxide`: an aqueous-CO₂ field, needed only when some PFT is carbon limited or an explicit
  calcifier (`nothing` otherwise)
- `nitrogen_to_carbon`, `zooplankton_phosphate_to_carbon`, `zooplankton_iron_to_carbon`: the fixed
  stoichiometric ratios
- `phosphate_quota_*`, `iron_quota_threshold_factor`, `*_silicon_quota`: variable growth-quota
  parameters
- `calcite_*`, `maximum_calcite_quota`: calcification parameters
- `loss_threshold_*_depth`, `zoo_loss_threshold_*_depth`, `aggregation_exponent`,
  `zoo_aggregation_exponent`: loss parameters
- `labile_fraction`, `nitrogen_to_DON_fraction`, `phosphorus_to_DOP_fraction`: routing of losses
  between the inorganic and dissolved-organic pools
- see the constructor definition for the full list of tunable parameters and their default
  values/units
"""
function ManyPhytoZoo(FT = Float64;
                      autotrophs  = MARBL_autotrophs(FT),
                      zooplankton = MARBL_zooplankton(FT),
                      carbon_dioxide = nothing,

                      zooplankton_phosphate_to_carbon = 1 / 117,
                      zooplankton_iron_to_carbon      = 3.0e-6,

                      nitrogen_to_carbon = 16 / 117,

                      growth_rate_regularisation   = 3.17e-8,   # 1/s
                      concentration_regularisation = 1.0e-8,    # mmol C/m³

                      labile_fraction            = 0.94,
                      nitrogen_to_DON_fraction   = 0.70,
                      phosphorus_to_DOP_fraction = 0.15,
                      aggregation_exponent       = 1.75,

                      phosphate_quota_slope       = 7.0,
                      phosphate_quota_intercept   = 5.571,
                      phosphate_to_carbon_maximum = 0.00854701,

                      iron_quota_threshold_factor = 10.0,

                      silicon_quota_reference        = 0.137,
                      maximum_silicon_quota          = 0.822,
                      minimum_silicon_quota          = 0.0457,
                      silicon_quota_threshold_factor = 6.0,

                      calcite_production_fraction   = 0.07,
                      maximum_calcifying_fraction   = 0.40,
                      calcite_bloom_threshold       = 2.5,     # mmol C/m³
                      calcite_temperature_threshold = 4.0,     # °C
                      calcite_temperature_minimum   = -2.0,    # °C
                      maximum_calcite_quota         = 0.4,

                      nitrogen_fixation_ratio = 1.25,

                      loss_threshold_full_depth = 80.0,        # m
                      loss_threshold_zero_depth = 120.0,       # m

                      zoo_aggregation_exponent      = 1.5,
                      zoo_loss_threshold_full_depth = 110.0,   # m
                      zoo_loss_threshold_zero_depth = 150.0,   # m

                      grazed_calcite_remineralisation_fraction = 0.33,
                      grazed_silica_remineralisation_fraction  = 0.50,
                      calcite_ballast_minimum      = 0.40,
                      small_phyto_poc_factor       = 0.13,     # 1/(mmol C)
                      grazed_small_phyto_poc_limit = 0.36,

                      alkalinity_phosphate_ratio = 0.0)  # 0 omits the PO₄–alkalinity term

    # define the tracer tendency methods for THIS PFT set's tracer names
    manifest_plankton!(autotrophs, zooplankton)

    # the union of the per-PFT limiting nutrients; each PFT applies its own subset, so this is carried
    # only to satisfy `AbstractPlankton{LN}`
    LN = (:nitrate, :ammonia, :phosphate, :iron, :silicate)

    return ManyPhytoZoo{LN, keys(autotrophs), keys(zooplankton), typeof(autotrophs),
                        typeof(zooplankton), typeof(carbon_dioxide), FT}(
        autotrophs,
        zooplankton,
        carbon_dioxide,
        convert(FT, zooplankton_phosphate_to_carbon),
        convert(FT, zooplankton_iron_to_carbon),
        convert(FT, nitrogen_to_carbon),
        convert(FT, growth_rate_regularisation),
        convert(FT, concentration_regularisation),
        convert(FT, labile_fraction),
        convert(FT, nitrogen_to_DON_fraction),
        convert(FT, phosphorus_to_DOP_fraction),
        convert(FT, aggregation_exponent),
        convert(FT, phosphate_quota_slope),
        convert(FT, phosphate_quota_intercept),
        convert(FT, phosphate_to_carbon_maximum),
        convert(FT, iron_quota_threshold_factor),
        convert(FT, silicon_quota_reference),
        convert(FT, maximum_silicon_quota),
        convert(FT, minimum_silicon_quota),
        convert(FT, silicon_quota_threshold_factor),
        convert(FT, calcite_production_fraction),
        convert(FT, maximum_calcifying_fraction),
        convert(FT, calcite_bloom_threshold),
        convert(FT, calcite_temperature_threshold),
        convert(FT, calcite_temperature_minimum),
        convert(FT, maximum_calcite_quota),
        convert(FT, nitrogen_fixation_ratio),
        convert(FT, loss_threshold_full_depth),
        convert(FT, loss_threshold_zero_depth),
        convert(FT, zoo_aggregation_exponent),
        convert(FT, zoo_loss_threshold_full_depth),
        convert(FT, zoo_loss_threshold_zero_depth),
        convert(FT, grazed_calcite_remineralisation_fraction),
        convert(FT, grazed_silica_remineralisation_fraction),
        convert(FT, calcite_ballast_minimum),
        convert(FT, small_phyto_poc_factor),
        convert(FT, grazed_small_phyto_poc_limit),
        convert(FT, alkalinity_phosphate_ratio))
end
