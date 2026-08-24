#####
##### The CESM2.1 (MARBL) plankton functional types
#####
##### Parameter values are CESM2.1, which differs from CESM2.0 in the iron and silicon growth quotas.
##### These are the building blocks the `MARBL` preset assembles.
#####

"""
    MARBL_small_phyto([FT = Float64;] kwargs...)

The CESM2.1 small phytoplankton (`sp`) — an implicit calcifier. Any keyword overrides the preset
value; see [`Autotroph`](@ref) for the full list.
"""
MARBL_small_phyto(FT = Float64; kwargs...) = Autotroph(FT;
    imp_calcifier = true,
    nitrate = 0.25, ammonia = 0.01, phosphate = 0.01, DOP = 0.3, iron = 3.0e-5, silicate = 0.0,
    photosynthesis_rate_reference_per_day = 5.0,
    maximum_chlorophyll_to_nitrogen = 2.5,
    temperature_threshold = -10.0,
    minimum_aggregation_rate_per_day = 0.01,
    loss_threshold_concentration = 0.01,
    cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 30.0e-6,
    minimum_iron_quota = 2.5e-6,
    kwargs...)

"""
    MARBL_diatoms([FT = Float64;] kwargs...)

The CESM2.1 diatoms (`diat`) — a silicifier. Any keyword overrides the preset value; see
[`Autotroph`](@ref) for the full list.
"""
MARBL_diatoms(FT = Float64; kwargs...) = Autotroph(FT;
    silicifier = true,
    nitrate = 0.5, ammonia = 0.05, phosphate = 0.05, DOP = 0.5, iron = 7.0e-5, silicate = 0.7,
    photosynthesis_rate_reference_per_day = 5.0,
    initial_PI_slope_per_day = 0.28,
    maximum_chlorophyll_to_nitrogen = 4.0,
    temperature_threshold = -10.0,
    minimum_aggregation_rate_per_day = 0.02,
    loss_threshold_concentration = 0.02,
    cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 30.0e-6,
    minimum_iron_quota = 2.5e-6,
    kwargs...)

"""
    MARBL_diazotrophs([FT = Float64;] kwargs...)

The CESM2.1 diazotrophs (`diaz`) — a nitrogen fixer. Any keyword overrides the preset value; see
[`Autotroph`](@ref) for the full list.
"""
MARBL_diazotrophs(FT = Float64; kwargs...) = Autotroph(FT;
    nfixer = true,
    nitrate = 2.0, ammonia = 0.2, phosphate = 0.015, DOP = 0.075, iron = 4.5e-5, silicate = 0.0,
    photosynthesis_rate_reference_per_day = 2.5,
    maximum_chlorophyll_to_nitrogen = 2.5,
    temperature_threshold = 15.0,
    minimum_aggregation_rate_per_day = 0.01,
    loss_threshold_concentration = 0.02,
    cold_loss_threshold_concentration = 0.001,
    growth_iron_quota_reference = 60.0e-6,
    minimum_iron_quota = 2.5e-6,
    kwargs...)

"""
    MARBL_autotrophs([FT = Float64])

The base CESM2.1 autotroph set: small phytoplankton, diatoms and diazotrophs.
"""
MARBL_autotrophs(FT = Float64) = (sp   = MARBL_small_phyto(FT),
                                  diat = MARBL_diatoms(FT),
                                  diaz = MARBL_diazotrophs(FT))

"""
    MARBL_zooplankton([FT = Float64])

The base CESM2.1 zooplankton: a single predator (`zoo`) grazing small phytoplankton, diatoms and
diazotrophs as three single-member prey classes.
"""
function MARBL_zooplankton(FT = Float64)
    grazing = (sp   = Grazing(FT; prey = :sp,
                              maximum_grazing_rate_per_day = 3.3,
                              fraction_to_zooplankton = 0.3,
                              fraction_to_particulate = 0.0,
                              detritus_fraction = 0.12),

               diat = Grazing(FT; prey = :diat,
                              maximum_grazing_rate_per_day = 3.15,
                              fraction_to_zooplankton = 0.25,
                              fraction_to_particulate = 0.39,
                              detritus_fraction = 0.24),

               diaz = Grazing(FT; prey = :diaz,
                              maximum_grazing_rate_per_day = 3.3,
                              fraction_to_zooplankton = 0.3,
                              fraction_to_particulate = 0.1,
                              detritus_fraction = 0.12))

    return (zoo = Heterotroph(FT; grazing), )
end

#####
##### The +cocco variant
#####
##### A four-PFT set in which coccolithophores are the sole calcifier and the only carbon-limited PFT.
##### The small phytoplankton do NOT calcify here — the same name carries different traits than in the
##### base set, which is why traits live per instance rather than in a name table — and the shared
##### parameters take their own values too.
#####

"""
    MARBL_cocco_small_phyto([FT = Float64;] kwargs...)

The small phytoplankton of the +cocco set. Unlike [`MARBL_small_phyto`](@ref) they do not calcify.
"""
MARBL_cocco_small_phyto(FT = Float64; kwargs...) = Autotroph(FT;
    nitrate = 0.2, ammonia = 0.01, phosphate = 0.005, DOP = 0.3, iron = 3.0e-5, silicate = 0.0,
    photosynthesis_rate_reference_per_day = 4.4,
    initial_PI_slope_per_day = 0.35,
    maximum_chlorophyll_to_nitrogen = 2.5,
    temperature_threshold = -10.0,
    minimum_aggregation_rate_per_day = 0.01,
    loss_threshold_concentration = 0.01,
    cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 35.0e-6,
    minimum_iron_quota = 3.0e-6,
    kwargs...)

"""
    MARBL_cocco_diatoms([FT = Float64;] kwargs...)

The diatoms of the +cocco set.
"""
MARBL_cocco_diatoms(FT = Float64; kwargs...) = Autotroph(FT;
    silicifier = true,
    nitrate = 0.5, ammonia = 0.05, phosphate = 0.05, DOP = 0.5, iron = 8.0e-5, silicate = 1.8,
    photosynthesis_rate_reference_per_day = 5.0,
    initial_PI_slope_per_day = 0.39,
    maximum_chlorophyll_to_nitrogen = 4.0,
    temperature_threshold = -10.0,
    minimum_aggregation_rate_per_day = 0.02,
    loss_threshold_concentration = 0.02,
    cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 35.0e-6,
    minimum_iron_quota = 3.0e-6,
    kwargs...)

"""
    MARBL_cocco_diazotrophs([FT = Float64;] kwargs...)

The diazotrophs of the +cocco set.
"""
MARBL_cocco_diazotrophs(FT = Float64; kwargs...) = Autotroph(FT;
    nfixer = true,
    nitrate = 2.0, ammonia = 0.2, phosphate = 0.015, DOP = 0.1, iron = 4.5e-5, silicate = 0.0,
    photosynthesis_rate_reference_per_day = 2.2,
    initial_PI_slope_per_day = 0.39,
    maximum_chlorophyll_to_nitrogen = 2.5,
    temperature_threshold = 15.0,
    minimum_aggregation_rate_per_day = 0.01,
    loss_threshold_concentration = 0.02,
    cold_loss_threshold_concentration = 0.001,
    growth_iron_quota_reference = 70.0e-6,
    minimum_iron_quota = 6.0e-6,
    kwargs...)

"""
    MARBL_coccolithophores([FT = Float64;] kwargs...)

Coccolithophores — an explicit calcifier which is carbon limited and uses the Eppley temperature form.
"""
MARBL_coccolithophores(FT = Float64; kwargs...) = Autotroph(FT;
    exp_calcifier = true, carbon_limited = true, temperature_function = :power,
    nitrate = 0.2, ammonia = 0.01, phosphate = 0.006, DOP = 0.25, iron = 3.15e-5, silicate = 0.0,
    photosynthesis_rate_reference_per_day = 4.7,
    initial_PI_slope_per_day = 0.28,
    maximum_chlorophyll_to_nitrogen = 3.5,
    temperature_threshold = 0.0,
    minimum_aggregation_rate_per_day = 0.01,
    loss_threshold_concentration = 0.01,
    cold_loss_threshold_concentration = 0.0,
    growth_iron_quota_reference = 35.0e-6,
    minimum_iron_quota = 3.0e-6,
    carbon_dioxide_half_saturation = 1.0,
    kwargs...)

"""
    MARBL_cocco_autotrophs([FT = Float64])

The four-PFT +cocco autotroph set.
"""
MARBL_cocco_autotrophs(FT = Float64) = (sp    = MARBL_cocco_small_phyto(FT),
                                        diat  = MARBL_cocco_diatoms(FT),
                                        diaz  = MARBL_cocco_diazotrophs(FT),
                                        cocco = MARBL_coccolithophores(FT))

"""
    MARBL_cocco_zooplankton([FT = Float64])

The +cocco zooplankton: one predator grazing all four autotrophs, with a sigmoidal response to the
small phytoplankton.
"""
function MARBL_cocco_zooplankton(FT = Float64)
    grazing = (sp    = Grazing(FT; prey = :sp,
                               maximum_grazing_rate_per_day = 3.6,
                               grazing_half_saturation = 0.54,
                               fraction_to_zooplankton = 0.3,
                               fraction_to_particulate = 0.0,
                               detritus_fraction = 0.12,
                               sigmoidal = true),

               diat  = Grazing(FT; prey = :diat,
                               maximum_grazing_rate_per_day = 3.41,
                               grazing_half_saturation = 0.72,
                               fraction_to_zooplankton = 0.25,
                               fraction_to_particulate = 0.38,
                               detritus_fraction = 0.24),

               diaz  = Grazing(FT; prey = :diaz,
                               maximum_grazing_rate_per_day = 3.35,
                               grazing_half_saturation = 0.6,
                               fraction_to_zooplankton = 0.3,
                               fraction_to_particulate = 0.1,
                               detritus_fraction = 0.12),

               cocco = Grazing(FT; prey = :cocco,
                               maximum_grazing_rate_per_day = 2.95,
                               grazing_half_saturation = 0.854,
                               fraction_to_zooplankton = 0.25,
                               fraction_to_particulate = 0.3,
                               detritus_fraction = 0.18))

    return (zoo = Heterotroph(FT; grazing), )
end

"""
    MARBL_cocco_plankton(grid, [FT = Float64;] kwargs...)

The four-PFT +cocco [`ManyPhytoZoo`](@ref). Needs the `grid` to allocate the aqueous-CO₂ field the
coccolithophore rate laws read.
"""
function MARBL_cocco_plankton(grid, FT = Float64;
                              autotrophs  = MARBL_cocco_autotrophs(FT),
                              zooplankton = MARBL_cocco_zooplankton(FT),
                              grazed_calcite_remineralisation_fraction = 0.70,
                              maximum_calcite_quota = 2.0,
                              kwargs...)

    return ManyPhytoZoo(FT;
                        autotrophs,
                        zooplankton,
                        carbon_dioxide = CenterField(grid),
                        grazed_calcite_remineralisation_fraction,
                        maximum_calcite_quota,
                        kwargs...)
end
