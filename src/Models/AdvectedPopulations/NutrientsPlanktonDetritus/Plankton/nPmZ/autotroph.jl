#####
##### Autotrophs
#####
struct Autotroph{TR, KN, FT}
            nutrient_half_saturations :: KN  # mmol/m³; MARBL kNO3/kNH4/kPO4/kDOP/kFe/kSiO3 (NamedTuple)

        photosynthesis_rate_reference :: FT  # 1/s; MARBL PCref
                     initial_PI_slope :: FT  # mmol C m²/(mg Chl W s); MARBL alphaPI
      maximum_chlorophyll_to_nitrogen :: FT  # mg Chl/mmol N; MARBL thetaN_max

              temperature_sensitivity :: FT  # ; MARBL Q_10
                temperature_reference :: FT  # °C; MARBL Tref
                temperature_threshold :: FT  # °C; MARBL temp_thres

                linear_mortality_rate :: FT  # 1/s; MARBL mort
             quadratic_mortality_rate :: FT  # 1/s / (mmol C/m³); MARBL mort2
             maximum_aggregation_rate :: FT  # 1/s; MARBL agg_rate_max
             minimum_aggregation_rate :: FT  # 1/s; MARBL agg_rate_min
            particulate_loss_fraction :: FT  # ; MARBL loss_poc

         loss_threshold_concentration :: FT  # mmol C/m³; MARBL loss_thres
    cold_loss_threshold_concentration :: FT  # mmol C/m³; MARBL loss_thres2

          growth_iron_quota_reference :: FT  # mmol Fe/mmol C; MARBL gQfe_0
                   minimum_iron_quota :: FT  # mmol Fe/mmol C; MARBL gQfe_min

       carbon_dioxide_half_saturation :: FT  # mmol/m³; MARBL kCO2
end

@inline traits(::Autotroph{TR}) where TR = TR

@inline carbon_name(s)       = Symbol(s, :C)
@inline chlorophyll_name(s)  = Symbol(s, :Chl)
@inline phosphorus_name(s)   = Symbol(s, :P)
@inline iron_name(s)         = Symbol(s, :Fe)
@inline silicon_name(s)      = Symbol(s, :Si)
@inline calcite_name(s)      = Symbol(s, :CaCO₃)

function autotroph_tracer_names(s, traits)
    names = [carbon_name(s), chlorophyll_name(s), phosphorus_name(s), iron_name(s)]
    traits.silicifier && push!(names, silicon_name(s))
    (traits.imp_calcifier || traits.exp_calcifier) && push!(names, calcite_name(s))
    return Tuple(names)
end

@inline calcifier(::Autotroph{traits}) where traits = traits.imp_calcifier | traits.exp_calcifier

const DEFAULT_AUTOTROPH_TRAITS = (silicifier           = false,   # silicifier       (⇒ <sname>Si tracer, VSi, gQsi)
                                  imp_calcifier        = false,   # imp_calcifier    (⇒ §9.1 CaCO3form)
                                  exp_calcifier        = false,   # exp_calcifier    (⇒ §9.2 picpoc)
                                  carbon_limited       = false,   # is_carbon_limited (⇒ §7 VCO2 term)
                                  nfixer               = false,   # Nfixer           (⇒ VN = 1, §10 Nfix)
                                  temperature_function = :q10)    # temp_func_form_opt (:q10 | :power)
