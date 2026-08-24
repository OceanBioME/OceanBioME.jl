Adapt.adapt_structure(to, a::Autotroph{TR, KN, FT}) where {TR, KN, FT} =
    Autotroph{TR, typeof(adapt(to, a.nutrient_half_saturations)), FT}(
        adapt(to, a.nutrient_half_saturations),
        a.photosynthesis_rate_reference,
        a.initial_PI_slope,
        a.maximum_chlorophyll_to_nitrogen,
        a.temperature_sensitivity,
        a.temperature_reference,
        a.temperature_threshold,
        a.linear_mortality_rate,
        a.quadratic_mortality_rate,
        a.maximum_aggregation_rate,
        a.minimum_aggregation_rate,
        a.particulate_loss_fraction,
        a.loss_threshold_concentration,
        a.cold_loss_threshold_concentration,
        a.growth_iron_quota_reference,
        a.minimum_iron_quota,
        a.carbon_dioxide_half_saturation)

Adapt.adapt_structure(to, g::Grazing) = g   # isbits

Adapt.adapt_structure(to, z::Heterotroph{GR, FT}) where {GR, FT} =
    Heterotroph{typeof(adapt(to, z.grazing)), FT}(
        adapt(to, z.grazing),
        z.linear_mortality_rate,
        z.quadratic_mortality_rate,
        z.loss_threshold_concentration,
        z.temperature_sensitivity,
        z.temperature_reference)

function Adapt.adapt_structure(to, p::ManyPhytoZoo{LN, P, Z}) where {LN, P, Z}
    autotrophs     = adapt(to, p.autotrophs)
    zooplankton    = adapt(to, p.zooplankton)
    carbon_dioxide = adapt(to, p.carbon_dioxide)

    return ManyPhytoZoo{LN, P, Z, typeof(autotrophs), typeof(zooplankton),
                        typeof(carbon_dioxide), typeof(p.nitrogen_to_carbon)}(
        autotrophs,
        zooplankton,
        carbon_dioxide,
        p.zooplankton_phosphate_to_carbon,
        p.zooplankton_iron_to_carbon,
        p.nitrogen_to_carbon,
        p.growth_rate_regularisation,
        p.concentration_regularisation,
        p.labile_fraction,
        p.nitrogen_to_DON_fraction,
        p.phosphorus_to_DOP_fraction,
        p.aggregation_exponent,
        p.phosphate_quota_slope,
        p.phosphate_quota_intercept,
        p.phosphate_to_carbon_maximum,
        p.iron_quota_threshold_factor,
        p.silicon_quota_reference,
        p.maximum_silicon_quota,
        p.minimum_silicon_quota,
        p.silicon_quota_threshold_factor,
        p.calcite_production_fraction,
        p.maximum_calcifying_fraction,
        p.calcite_bloom_threshold,
        p.calcite_temperature_threshold,
        p.calcite_temperature_minimum,
        p.maximum_calcite_quota,
        p.nitrogen_fixation_ratio,
        p.loss_threshold_full_depth,
        p.loss_threshold_zero_depth,
        p.zoo_aggregation_exponent,
        p.zoo_loss_threshold_full_depth,
        p.zoo_loss_threshold_zero_depth,
        p.grazed_calcite_remineralisation_fraction,
        p.grazed_silica_remineralisation_fraction,
        p.calcite_ballast_minimum,
        p.small_phyto_poc_factor,
        p.grazed_small_phyto_poc_limit,
        p.alkalinity_phosphate_ratio)
end

Base.summary(::ManyPhytoZoo{LN, P, Z}) where {LN, P, Z} =
    "ManyPhytoZoo ($(length(P)) autotroph$(length(P) == 1 ? "" : "s") + $(length(Z)) zooplankton)"

# describe each PFT by the traits it actually carries, so a user's PFT set reads as well as a standard one
function autotroph_description(p::ManyPhytoZoo, s)
    t = traits(p, Val(s))

    roles = String[]

    t.imp_calcifier  && push!(roles, "implicit calcifier")
    t.exp_calcifier  && push!(roles, "explicit calcifier")
    t.silicifier     && push!(roles, "silicifier")
    t.nfixer         && push!(roles, "N-fixer")
    t.carbon_limited && push!(roles, "carbon limited")

    return "$s: " * (isempty(roles) ? "generic autotroph" : join(roles, ", "))
end

# "zoo eats (sp), (diat, cocco)" — one entry per prey class of that predator
function heterotroph_description(p::ManyPhytoZoo, z)
    classes = values(getproperty(p.zooplankton, z).grazing)

    return "$z eats " * join((string(prey_snames(g)) for g in classes), ", ")
end

function Base.show(io::IO, p::ManyPhytoZoo{LN, P, Z}) where {LN, P, Z}
    msg  = summary(p) * "\n"
    msg *= "├── autotrophs: " * join((autotroph_description(p, s) for s in P), " │ ") * "\n"
    msg *= "├── zooplankton: " * join((heterotroph_description(p, z) for z in Z), " │ ") * "\n"
    msg *= "├── N:C (fixed) = $(p.nitrogen_to_carbon)\n"
    msg *= "├── zooplankton P:C = $(p.zooplankton_phosphate_to_carbon), Fe:C = $(p.zooplankton_iron_to_carbon)\n"
    msg *= "└── tracers: $(required_biogeochemical_tracers(p))"

    print(io, msg)

    return nothing
end
