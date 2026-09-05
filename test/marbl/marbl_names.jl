#####
##### The names the MARBL Fortran-comparison harnesses use, aliased onto the OceanBioME components.
#####
##### The harnesses were written against the standalone prototype, whose types carried MARBL-specific
##### names. The components are now general OceanBioME ones with MARBL as a parameter set, so this file
##### is the single place the two vocabularies meet — everything else in `test/marbl/` reads as the
##### harness always did.
#####

using OceanBioME
using Oceananigans.Fields: CenterField

const NPDM = OceanBioME.Models.NutrientsPlanktonDetritusModels
const MPK  = NPDM.PlanktonModels.MARBLPlankton
const DET  = NPDM.DetritusModels.MultiElement
const NUT  = NPDM.NutrientsModels
const OXY  = NPDM.OxygenModels
const ICM  = NPDM.InorganicCarbonModels
const SED  = OceanBioME.Models.SedimentModels

# components
using .NPDM.PlanktonModels.MARBLPlankton: ManyPhytoZoo as MARBLPlankton,
                                          Autotroph as MARBLAutotroph,
                                          Heterotroph as MARBLZooplankton,
                                          Grazing as MARBLGrazing
using .NPDM.DetritusModels: MultiElementRefractoryDissolvedParticulate as MARBLDetritus,
                            MultiElementRefractoryDissolved as MARBLBallastDetritus
using .NPDM.DetritusModels.MultiElement: Ballast as MARBLBallast
using .NPDM.OxygenModels: RedoxOxygen as MARBLOxygen
using .NPDM.InorganicCarbonModels: ImplicitExplicitCalcite
using .NPDM.NutrientsModels: ComplexedIron
using .OceanBioME.Models.SedimentModels: BurialDenitrification as MARBLSediment

MARBLPlankton_cocco(grid, FT = Float64; kwargs...) = MPK.MARBL_cocco_plankton(grid, FT; kwargs...)
MARBLSedimentModel(grid; kwargs...) = SED.BurialDenitrificationSediment(grid; kwargs...)

# the CESM2.1 parameter sets, under the names the harnesses use
const default_small_phytoplankton = MPK.MARBL_small_phyto
const default_diatoms             = MPK.MARBL_diatoms
const default_diazotrophs         = MPK.MARBL_diazotrophs
const default_coccolithophores    = MPK.MARBL_coccolithophores
const cocco_small_phytoplankton   = MPK.MARBL_cocco_small_phyto
const cocco_diatoms               = MPK.MARBL_cocco_diatoms
const cocco_diazotrophs           = MPK.MARBL_cocco_diazotrophs

# rate functions and helpers reached for directly
for (mod, names) in ((MPK, (:traits, :carbon_name, :chlorophyll_name, :phosphorus_name, :iron_name,
                            :silicon_name, :calcite_name, :autotroph_tracer_names,
                            :temperature_scaling, :max_specific_rate, :light_limited_rate,
                            :photoacclimation_rate, :growth_phosphorus_quota, :growth_iron_quota,
                            :growth_silicon_quota, :aggregation_rate, :nitrogen_fixation,
                            :grazing_rate, :calcifier_graze_poc_fraction, :zoo_loss_rate,
                            :calcite_formation, :photosynthesis,
                            :autotroph_aggregation, :autotroph_linear_loss, :loss_to_poc,
                            :graze_to_doc, :graze_to_poc, :graze_to_predator, :graze_to_zooplankton,
                            :grazing_loss, :zoo_loss_to_doc, :zoo_loss_to_poc,
                            :zooplankton_grazing_gain, :zooplankton_loss)),
                     (NUT, (:iron_speciation, :iron_scavenging_flux, :ligand_scavenging_flux,
                            :sinking_mass)),
                     (OXY, (:water_column_nitrification, :water_column_denitrification)),
                     (SED, (:organic_carbon_burial_fraction, :organic_nitrogen_burial_fraction,
                            :calcite_burial_fraction, :iron_burial_fraction,
                            :sedimentary_denitrification, :other_remineralisation,
                            :denitrification_response, :coupled_tracers)),
                     (ICM, (:particulate_calcium_carbonate_production, :carbonate_concentration)))
    for n in names
        isdefined(mod, n) && @eval const $n = $mod.$n
    end
end

#####
##### Signatures that changed in the port
#####

# the water-column and benthic oxygen consumption were one function in the prototype and are two here
function oxygen_consumption end

@inline oxygen_consumption(i, j, k, grid, bgc, fields, aux) =
    OXY.oxygen_consumption(i, j, k, grid, bgc, fields, aux)

@inline oxygen_consumption(s::SED.BurialDenitrification, args...) = SED.oxygen_consumption(s, args...)

# these gained an argument: the plankton, and the iron cycle whose parameters replaced hardcoded constants
@inline oxygen_production_total(i, j, k, grid, bgc, fields, aux) =
    MPK.oxygen_production_total(i, j, k, grid, bgc.plankton, bgc, fields, aux)

@inline ligand_photochemistry_rate(i, j, k, grid, aux) =
    NUT.ligand_photochemistry_rate(i, j, k, grid, ComplexedIron(), aux)

#####
##### Tracer names
#####

# our tracer names carry subscripts; MARBL's baseline variables are ASCII
ascii_name(n) = replace(String(n), "₀"=>"0", "₁"=>"1", "₂"=>"2", "₃"=>"3", "₄"=>"4")

# and the inverse, for turning a MARBL variable name back into one of our tracer symbols. Calcite is the
# only plankton tracer suffix that carries a subscript; the rest (C, Chl, P, Fe, Si) are already plain.
tracer_symbol(n) = Symbol(replace(String(n), "CaCO3" => "CaCO₃"))

# temperature is a physics tracer, not a biogeochemical one, so a harness must set it itself
const PHYSICS_TRACERS = (:T, )
