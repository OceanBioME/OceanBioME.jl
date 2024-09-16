"""
Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model.

Tracers
=======

* Nano-phytoplankton: P (μmolC/L)
* Diatoms: D (μmolC/L)
* Zooplankton: Z (μmolC/L)
* Mesozooplankton: M (μmolC/L)
* Chlorophyll in nano-phytoplankton: Pᶜʰˡ (μgChl/L)
* Chlorophyll in diatoms: Dᶜʰˡ (μgChl/L)
* Iron in nano-phytoplanktons: Pᶠᵉ (nmolFe/L) 
* Iron in diatoms: Dᶠᵉ (nmolFe/L) 
* Silicon in diatoms: Dˢⁱ (μmolSi/L)

* Dissolved organic carbon: DOC (μmolC/L)
* Small sinking particles : POC (μmolC/L)
* Large sinking particles: GOC (μmolC/L)
* Iron in small particles: SFe (nmolFe/L) 
* Iron in large particles: BFe (nmolFe/L) 
* Silicate in large particles : PSi (μmolSi/L)
* Nitrates: NO₃ (μmolN/L)
* Ammonium: NH₄ (μmolN/L)
* Phosphate: PO₄ (μmolP/L)
* Dissolved iron: Fe (nmolFe/L)
* Silicate: Si (μmolSi/L)
* Calcite: CaCO₃ (μmolC/L)
* Dissolved oxygen: O₂ (μmolO₂/L)
* Dissolved inorganic carbon: DIC (μmolC/L)
* Total alkalinity: Alk (μmolN/L)


Required submodels
==================
# you will need something like this, they use a different PAR model but I wouldn't worry about that for now, you might also need temperatre and salinity (not sure)
* Photosynthetically available radiation: PAR (W/m²)

"""
module PISCESModel

export PISCES

using Oceananigans.Units

using Oceananigans: KernelFunctionOperation
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField, ConstantField, Center, Face

using OceanBioME.Light: MultiBandPhotosyntheticallyActiveRadiation, default_surface_PAR, compute_euphotic_depth!
using OceanBioME: setup_velocity_fields, show_sinking_velocities, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.BoxModels: BoxModel
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: set!
using Oceananigans.Grids: φnodes, RectilinearGrid

import OceanBioME: redfield, conserved_tracers, maximum_sinking_velocity, chlorophyll

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity,
                                     biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt
import Base: show, summary

import OceanBioME.Models.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers
struct PISCES{NP, DP, SZ, BZ, DM, PM, NI, FE, SI, OX, PO, CA, CE, FT, LA, DL, ML, EU, MS, DD, VD, MP, CC, CS, SS} <: AbstractContinuousFormBiogeochemistry
                        nanophytoplankton :: NP
                                  diatoms :: DP

                         microzooplankton :: SZ
                          mesozooplankton :: BZ

                 dissolved_organic_matter :: DM
               particulate_organic_matter :: PM

                                 nitrogen :: NI
                                     iron :: FE
                                 silicate :: SI
                                   oxygen :: OX
                                phosphate :: PO

                                  calcite :: CA
                            carbon_system :: CE

                   first_anoxia_threshold :: FT
                  second_anoxia_threshold :: FT

                  nitrogen_redfield_ratio :: FT
                 phosphate_redfield_ratio :: FT
                      iron_redfield_ratio :: FT

                        mixed_layer_shear :: FT
                         background_shear :: FT

                                 latitude :: LA
                               day_length :: DL

                        mixed_layer_depth :: ML
                           euphotic_depth :: EU
                  yearly_maximum_silicate :: MS
                          dust_deposition :: DD

    mean_mixed_layer_vertical_diffusivity :: VD
                   mean_mixed_layer_light :: MP

                         carbon_chemistry :: CC
                       calcite_saturation :: CS

                       sinking_velocities :: SS
end

const NANO_PHYTO = Union{Val{:P}, Val{:PChl}, Val{:PFe}}
const DIATOMS    = Union{Val{:D}, Val{:DChl}, Val{:DFe}, Val{:DSi}}
const PARTICLES = Union{Val{:POC}, Val{:SFe}, Val{:GOC}, Val{:BFe}, Val{:PSi}}
const NITROGEN = Union{Val{:NO₃}, Val{:NH₄}}
const CARBON_SYSTEM = Union{Val{:DIC}, Val{:Alk}}

(bgc::PISCES)(val_name::NANO_PHYTO, args...)    = bgc.nanophytoplankton(val_name, bgc, args...)
(bgc::PISCES)(val_name::DIATOMS, args...)       = bgc.diatoms(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:Z}, args...)       = bgc.microzooplankton(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:M}, args...)       = bgc.mesozooplankton(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:DOC}, args...)     = bgc.dissolved_organic_matter(val_name, bgc, args...)
(bgc::PISCES)(val_name::PARTICLES, args...)     = bgc.particulate_organic_matter(val_name, bgc, args...)
(bgc::PISCES)(val_name::NITROGEN, args...)      = bgc.nitrogen(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:Fe}, args...)      = bgc.iron(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:Si}, args...)      = bgc.silicate(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:CaCO₃}, args...)   = bgc.calcite(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:O₂}, args...)      = bgc.oxygen(val_name, bgc, args...)
(bgc::PISCES)(val_name::Val{:PO₄}, args...)     = bgc.phosphate(val_name, bgc, args...)
(bgc::PISCES)(val_name::CARBON_SYSTEM, args...) = bgc.carbon_system(val_name, bgc, args...)


@inline biogeochemical_auxiliary_fields(bgc::PISCES) = 
    (zₘₓₗ = bgc.mixed_layer_depth, 
     zₑᵤ = bgc.euphotic_depth, 
     Si̅ = bgc.yearly_maximum_silicate, 
     D_dust = bgc.dust_deposition, 
     Ω = bgc.calcite_saturation,
     κ = bgc.mean_mixed_layer_vertical_diffusivity,
     mixed_layer_PAR = bgc.mean_mixed_layer_light)

@inline required_biogeochemical_tracers(::PISCES) = 
    (:P, :D, :Z, :M, :PChl, :DChl, :PFe, :DFe, :DSi, 
     :DOC, :POC, :GOC, :SFe, :BFe, :PSi, # its really silly that this is called PSi when DSi also exists
     :NO₃, :NH₄, :PO₄, :Fe, :Si, 
     :CaCO₃, :DIC, :Alk, :O₂, :T, :S)

@inline required_biogeochemical_auxiliary_fields(::PISCES) =
    (:zₘₓₗ, :zₑᵤ, :Si̅, :D_dust, :Ω, :κ, :mixed_layer_PAR, :PAR, :PAR₁, :PAR₂, :PAR₃)

const small_particle_components = Union{Val{:POC}, Val{:SFe}}
const large_particle_components = Union{Val{:GOC}, Val{:BFe}, Val{:PSi}, Val{:CaCO₃}} 

biogeochemical_drift_velocity(bgc::PISCES, ::small_particle_components) = (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.POC)
biogeochemical_drift_velocity(bgc::PISCES, ::large_particle_components) = (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities.GOC)

summary(::PISCES{FT}) where {FT} = string("PISCES{$FT}") 

show(io::IO, model::PISCES) = print(io, string("Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model")) # maybe add some more info here

include("common.jl")
include("phytoplankton.jl")
include("zooplankton.jl")
include("dissolved_organic_matter.jl")
include("particulate_organic_matter.jl")
include("nitrate_ammonia.jl")
include("phosphates.jl")
include("iron.jl")
include("silicon.jl")
include("calcite.jl")
include("carbonate_system.jl")
include("oxygen.jl")
include("mean_mixed_layer_properties.jl")
include("compute_calcite_saturation.jl")
include("update_state.jl")
include("coupling_utils.jl")

# to change to new production change `NutrientLimitedProduction` for `GrowthRespirationLimitedProduction`
function PISCES(; grid,
                  nanophytoplankton = 
                    Phytoplankton(growth_rate = NutrientLimitedProduction(dark_tollerance = 3days),
                                  nutrient_limitation = 
                                    NitrogenIronPhosphateSilicateLimitation(minimum_ammonium_half_saturation = 0.013,
                                                                            minimum_nitrate_half_saturation = 0.13, 
                                                                            minimum_phosphate_half_saturation = 0.8,
                                                                            half_saturation_for_iron_uptake = 1.0,
                                                                            silicate_limited = false),
                                  blue_light_absorption = 2.1, 
                                  green_light_absorption = 0.42, 
                                  red_light_absorption = 0.4,
                                  maximum_quadratic_mortality = 0.0,
                                  maximum_chlorophyll_ratio = 0.033),

                  diatoms = 
                    Phytoplankton(growth_rate = NutrientLimitedProduction(dark_tollerance = 4days),
                                  nutrient_limitation = 
                                    NitrogenIronPhosphateSilicateLimitation(minimum_ammonium_half_saturation = 0.039,
                                                                            minimum_nitrate_half_saturation = 0.39, 
                                                                            minimum_phosphate_half_saturation = 2.4,
                                                                            half_saturation_for_iron_uptake = 3.0,
                                                                            silicate_limited = true),
                                  blue_light_absorption = 1.6, 
                                  green_light_absorption = 0.69, 
                                  red_light_absorption = 0.7,
                                  maximum_quadratic_mortality = 0.03/day,
                                  maximum_chlorophyll_ratio = 0.05),
                  
                  microzooplankton = Zooplankton(maximum_grazing_rate = 3/day,
                                                 preference_for_nanophytoplankton = 1.0,
                                                 preference_for_diatoms = 0.5,
                                                 preference_for_particulates = 0.1,
                                                 preference_for_zooplankton = 0.0,
                                                 quadratic_mortality = 0.004/day,
                                                 linear_mortality = 0.03/day,
                                                 minimum_growth_efficiency = 0.3,
                                                 maximum_flux_feeding_rate = 0.0,
                                                 undissolved_calcite_fraction = 0.5),

                  mesozooplankton = Zooplankton(maximum_grazing_rate = 0.75/day,
                                                preference_for_nanophytoplankton = 0.3,
                                                preference_for_diatoms = 1.0,
                                                preference_for_particulates = 0.3,
                                                preference_for_zooplankton = 1.0,
                                                quadratic_mortality = 0.03/day,
                                                linear_mortality = 0.005/day,
                                                minimum_growth_efficiency = 0.35,
                                                # not documented but the below must implicitly contain a factor of second/day
                                                # to be consistent in the NEMO namelist to go from this * mol / L * m/s to mol / L / day
                                                maximum_flux_feeding_rate = 2e3 / 1e6 / day, # (day * meter/s * mol/L)^-1 to (meter * μ mol/L)^-1
                                                undissolved_calcite_fraction = 0.75),
                  
                  dissolved_organic_matter = DissolvedOrganicMatter(),
                  particulate_organic_matter = TwoCompartementParticulateOrganicMatter(),
                  
                  nitrogen = NitrateAmmonia(),
                  iron = SimpleIron(),
                  silicate = Silicate(),
                  oxygen = Oxygen(),
                  phosphate = Phosphate(),
                  
                  calcite = Calcite(),
                  carbon_system = CarbonateSystem(),

                  # from Aumount 2005 rather than 2015 since it doesn't work the other way around
                  first_anoxia_thresehold = 6.0,
                  second_anoxia_thresehold = 1.0,

                  nitrogen_redfield_ratio = 16/122,
                  phosphate_redfield_ratio = 1/122,
                  iron_redfield_ratio = 10^-3,
                  
                  mixed_layer_shear = 1.0/day,
                  background_shear = 0.01/day, 
                  
                  latitude = PrescribedLatitude(45),
                  day_length = day_length_function,
                  
                  mixed_layer_depth = Field{Center, Center, Nothing}(grid),
                  euphotic_depth = Field{Center, Center, Nothing}(grid),

                  yearly_maximum_silicate = ConstantField(7.5),
                  dust_deposition = ZeroField(),

                  mean_mixed_layer_vertical_diffusivity = Field{Center, Center, Nothing}(grid),
                  mean_mixed_layer_light = Field{Center, Center, Nothing}(grid),

                  carbon_chemistry = CarbonChemistry(),
                  calcite_saturation = CenterField(grid),

                  surface_photosynthetically_active_radiation = default_surface_PAR,

                  light_attenuation =
                    MultiBandPhotosyntheticallyActiveRadiation(; grid, 
                                                                 surface_PAR = surface_photosynthetically_active_radiation),

                  sinking_speeds = (POC = 2/day, 
                                    # might be more efficient to just precompute this
                                    GOC = KernelFunctionOperation{Center, Center, Face}(DepthDependantSinkingSpeed(), 
                                                                                        grid, 
                                                                                        mixed_layer_depth, 
                                                                                        euphotic_depth)),
                  open_bottom = true,

                  scale_negatives = false,
                  invalid_fill_value = NaN,
                  
                  sediment = nothing,
                  particles = nothing,
                  modifiers = nothing)

    @warn "This implementation of PISCES is in early development and has not yet been fully validated"

    if !isnothing(sediment) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    sinking_velocities = merge(sinking_velocities, (; grid)) # we need to interpolate the fields so we need this for flux feeding inside a kernel - this might cause problems...

    if (latitude isa PrescribedLatitude) & !(grid isa RectilinearGrid)
        φ = φnodes(grid, Center(), Center(), Center())

        @warn "A latitude of $latitude was given but the grid has its own latitude ($(minimum(φ)), $(maximum(φ))) so the prescribed value is ignored"

        latitude = nothing 
    elseif (latitude isa ModelLatitude) & (grid isa RectilinearGrid)
        throw(ArgumentError("You must prescribe a latitude when using a `RectilinearGrid`"))
    end

    # just incase we're in the default state with no closure model
    # this highlights that the darkness term for phytoplankton growth is obviously wrong because not all phytoplankon
    # cells spend an infinite amount of time in the dark if the diffusivity is zero, it should depend on where they are...
    if !(mean_mixed_layer_vertical_diffusivity isa ConstantField)
        set!(mean_mixed_layer_vertical_diffusivity, 1)
    end

    underlying_biogeochemistry = PISCES(nanophytoplankton, diatoms,
                                        microzooplankton, mesozooplankton,
                                        dissolved_organic_matter, particulate_organic_matter,
                                        nitrogen, iron, silicate, oxygen, phosphate,
                                        calcite, carbon_system, 
                                        first_anoxia_thresehold, second_anoxia_thresehold,
                                        nitrogen_redfield_ratio, phosphate_redfield_ratio, iron_redfield_ratio,
                                        mixed_layer_shear, background_shear,
                                        latitude, day_length,
                                        mixed_layer_depth, euphotic_depth,
                                        yearly_maximum_silicate, dust_deposition,
                                        mean_mixed_layer_vertical_diffusivity, 
                                        mean_mixed_layer_light,
                                        carbon_chemistry, calcite_saturation,
                                        sinking_velocities)

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


function total_carbon_tendency(bgc, x, y, z, t,
                               P, D, Z, M, 
                               PChl, DChl, PFe, DFe, DSi,
                               DOC, POC, GOC, 
                               SFe, BFe, PSi, 
                               NO₃, NH₄, PO₄, Fe, Si, 
                               CaCO₃, DIC, Alk, 
                               O₂, T, S,
                               zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    # nano phyto
    phyto = bgc.nanophytoplankton
    val_name = Val(:P)

    # mortality
    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, I, zₘₓₗ, L)

    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.microzooplankton, P, D, Z, POC, T)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = gZ * Z + gM * M

    @info production, linear_mortality, quadratic_mortality, grazing#production - linear_mortality - quadratic_mortality - grazing

    # diatoms 
    phyto = bgc.diatoms
    val_name = Val(:D)

    # mortality
    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, I, zₘₓₗ, L)

    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.microzooplankton, P, D, Z, POC, T)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = gZ * Z + gM * M

    @info production, linear_mortality, quadratic_mortality, grazing#production - linear_mortality - quadratic_mortality - grazing

    # micro zoo
    zoo = bgc.microzooplankton
    val_name = Val(:Z)

    I = zooplankton_concentration(val_name, Z, M)

    # grazing
    total_specific_grazing, gP, gD, gPOC, gZ = specific_grazing(zoo, P, D, Z, POC, T)

    grazing = total_specific_grazing * I

    # flux feeding
    grid = bgc.sinking_velocities.grid

    small_flux_feeding = specific_flux_feeding(zoo, x, y, z, POC, T, bgc.sinking_velocities.POC, grid)
    large_flux_feeding = specific_flux_feeding(zoo, x, y, z, GOC, T, bgc.sinking_velocities.GOC, grid)

    flux_feeding = (small_flux_feeding + large_flux_feeding) * I

    # grazing mortality
    specific_grazing_mortality = zooplankton_grazing_mortality(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing_mortality = specific_grazing_mortality * M

    # mortality
    total_mortality = mortality(zoo, bgc, I, O₂, T)

    growth_efficiency = grazing_growth_efficiency(zoo, P, D, PFe, DFe, POC, SFe, total_specific_grazing, gP, gD, gPOC, gZ)

    @info  growth_efficiency .* (grazing, flux_feeding), grazing_mortality, total_mortality#growth_efficiency * (grazing + flux_feeding) - grazing_mortality - total_mortality

    # mesozoo
    # micro zoo
    zoo = bgc.microzooplankton
    val_name = Val(:Z)

    I = zooplankton_concentration(val_name, Z, M)

    # grazing
    total_specific_grazing, gP, gD, gPOC, gZ = specific_grazing(zoo, P, D, Z, POC, T)

    grazing = total_specific_grazing * I

    # flux feeding
    grid = bgc.sinking_velocities.grid

    small_flux_feeding = specific_flux_feeding(zoo, x, y, z, POC, T, bgc.sinking_velocities.POC, grid)
    large_flux_feeding = specific_flux_feeding(zoo, x, y, z, GOC, T, bgc.sinking_velocities.GOC, grid)

    flux_feeding = (small_flux_feeding + large_flux_feeding) * I

    # grazing mortality
    specific_grazing_mortality = zooplankton_grazing_mortality(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing_mortality = specific_grazing_mortality * M

    # mortality
    total_mortality = mortality(zoo, bgc, I, O₂, T)

    growth_efficiency = grazing_growth_efficiency(zoo, P, D, PFe, DFe, POC, SFe, total_specific_grazing, gP, gD, gPOC, gZ)

    @info  growth_efficiency .* (grazing, flux_feeding), grazing_mortality, total_mortality#growth_efficiency * (grazing + flux_feeding) - grazing_mortality - total_mortality

    # DOC
    dom = bgc.dissolved_organic_matter
    val_name = Val(DOC)

    nanophytoplankton_exudation = dissolved_exudate(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    diatom_exudation = dissolved_exudate(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    phytoplankton_exudation = nanophytoplankton_exudation + diatom_exudation

    particulate_degredation = specific_degredation_rate(bgc.particulate_organic_matter, bgc, O₂, T) * POC

    respiration_product = dissolved_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    microzooplankton_grazing_waste = specific_dissolved_grazing_waste(bgc.microzooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * Z
    mesozooplankton_grazing_waste  = specific_dissolved_grazing_waste(bgc.mesozooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * M

    grazing_waste = microzooplankton_grazing_waste + mesozooplankton_grazing_waste

    degredation = bacterial_degradation(dom, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    aggregation_to_particles, = aggregation(dom, bgc, z, DOC, POC, GOC, zₘₓₗ)

    @info phytoplankton_exudation, particulate_degredation, respiration_product, grazing_waste, degredation, aggregation_to_particles# phytoplankton_exudation + particulate_degredation + respiration_product + grazing_waste - degredation - aggregation_to_particles

    # poc
    poc = bgc.particulate_organic_matter
    val_name = Val(:POC)
    grazing_waste = specific_non_assimilated_waste(bgc.microzooplankton, bgc, x, y, z, P, D, Z, POC, GOC, T) * Z

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)

    nanophytoplankton_mortality = (1 - 0.5 * R_CaCO₃) * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality)

    diatom_linear_mortality, = mortality(bgc.diatoms, bgc, z, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)

    diatom_mortality = 0.5 * diatom_linear_mortality

    microzooplankton_mortality = mortality(bgc.microzooplankton, bgc, Z, O₂, T)

    # degredation
    λ = specific_degredation_rate(poc, bgc, O₂, T)

    large_particle_degredation = λ * GOC
    degredation = λ * POC

    # grazing
    _, _, _, microzooplankton_grazing = specific_grazing(bgc.microzooplankton, P, D, Z, POC, T)
    _, _, _, mesozooplankton_grazing = specific_grazing(bgc.mesozooplankton, P, D, Z, POC, T)

    grid = bgc.sinking_velocities.grid

    small_flux_feeding = specific_flux_feeding(bgc.mesozooplankton, x, y, z, POC, T, bgc.sinking_velocities.POC, grid)

    grazing = microzooplankton_grazing * Z + (mesozooplankton_grazing + small_flux_feeding) * M

    # aggregation
    _, Φ₁, _, Φ₃ = aggregation(bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, zₘₓₗ)
    dissolved_aggregation = Φ₁ + Φ₃
    
    aggregation_to_large = aggregation(poc, bgc, z, POC, GOC, zₘₓₗ)

    total_aggregation = dissolved_aggregation - aggregation_to_large

    @info (grazing_waste 
            ,nanophytoplankton_mortality ,diatom_mortality ,microzooplankton_mortality 
            ,large_particle_degredation ,total_aggregation 
        , grazing, degredation)

    # goc
    poc = bgc.particulate_organic_matter
    val_name = Val(:GOC)
    grazing_waste = specific_non_assimilated_waste(bgc.mesozooplankton, bgc, x, y, z, P, D, Z, POC, GOC, T) * M

    # mortality terms
    R_CaCO₃ = rain_ratio(bgc.calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    nanophytoplankton_linear_mortality, nanophytoplankton_quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)

    nanophytoplankton_mortality = 0.5 * R_CaCO₃ * (nanophytoplankton_linear_mortality + nanophytoplankton_quadratic_mortality)

    diatom_linear_mortality, diatom_quadratic_mortality = mortality(bgc.diatoms, bgc, z, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)

    diatom_mortality = 0.5 * diatom_linear_mortality + diatom_quadratic_mortality

    mesozooplankton_mortality = mortality(bgc.mesozooplankton, bgc, M, O₂, T)

    # degredation
    λ = specific_degredation_rate(poc, bgc, O₂, T)

    degredation = λ * GOC

    # grazing
    grid = bgc.sinking_velocities.grid
    grazing = specific_flux_feeding(bgc.mesozooplankton, x, y, z, GOC, T, bgc.sinking_velocities.GOC, grid) * M

    # aggregation
    _, _, dissolved_aggregation = aggregation(bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, zₘₓₗ)
    
    small_particle_aggregation = aggregation(poc, bgc, z, POC, GOC, zₘₓₗ)

    total_aggregation = dissolved_aggregation + small_particle_aggregation

    # fecal pelet prodiction
    fecal_pelet_production = upper_trophic_fecal_product(bgc.mesozooplankton, M, T)

    #@info grazing_waste, nanophytoplankton_mortality , diatom_mortality , mesozooplankton_mortality, total_aggregation , fecal_pelet_production, grazing , degredation
    @info (grazing_waste 
    ,nanophytoplankton_mortality ,diatom_mortality ,mesozooplankton_mortality 
    ,total_aggregation ,fecal_pelet_production
    ,grazing 
    ,degredation)
    #=return (grazing_waste 
            + nanophytoplankton_mortality + diatom_mortality + mesozooplankton_mortality 
            + total_aggregation + fecal_pelet_production
            - grazing 
            - degredation)=#

    # DIC

    microzooplankton_respiration = specific_inorganic_grazing_waste(bgc.microzooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * Z
    mesozooplankton_respiration  = specific_inorganic_grazing_waste(bgc.mesozooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, T) * M

    zooplankton_respiration = microzooplankton_respiration + mesozooplankton_respiration

    upper_trophic_respiration = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    dissolved_degredation = bacterial_degradation(bgc.dissolved_organic_matter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    calcite_diss = calcite_dissolution(bgc.calcite, CaCO₃, Ω)

    calcite_prod = calcite_production(bgc.calcite, bgc, z, P, D, PChl, PFe, Z, M, POC, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    calcite = calcite_diss - calcite_prod

    nanophytoplankton_consumption = total_production(bgc.nanophytoplankton, bgc, y, t, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    diatom_consumption            = total_production(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    consumption = nanophytoplankton_consumption + diatom_consumption
    
    @info zooplankton_respiration, upper_trophic_respiration, dissolved_degredation, calcite, consumption#zooplankton_respiration + upper_trophic_respiration + dissolved_degredation + calcite - consumption

    # CaCO\_3   
    calcite = bgc.calcite
    production = calcite_production(calcite, bgc, z, P, D, PChl, PFe, Z, M, POC, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    dissolution = calcite_dissolution(calcite, CaCO₃, Ω)
    
    @info production, dissolution, production - dissolution#production - dissolution
end


function zooplankton_mortality_terms(bgc, x, y, z, t,
                               P, D, Z, M, 
                               PChl, DChl, PFe, DFe, DSi,
                               DOC, POC, GOC, 
                               SFe, BFe, PSi, 
                               NO₃, NH₄, PO₄, Fe, Si, 
                               CaCO₃, DIC, Alk, 
                               O₂, T, S,
                               zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    val_name = Val(:Z)
    zoo = bgc.microzooplankton
    I = zooplankton_concentration(val_name, Z, M)

    # mortality
    Z_mortality = mortality(zoo, bgc, I, O₂, T)

    val_name = Val(:M)
    zoo = bgc.mesozooplankton
    I = zooplankton_concentration(val_name, Z, M)

    # mortality
    M_mortality = mortality(zoo, bgc, I, O₂, T)

    # DOC

    DOC_respiration_product = dissolved_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    # DIC
    DIC_respiration_product = inorganic_upper_trophic_respiration_product(bgc.mesozooplankton, M, T)

    # POC
    POC_microzooplankton_mortality = mortality(bgc.microzooplankton, bgc, Z, O₂, T)

    # GOC
    GOC_mesozooplankton_mortality = linear_mortality(bgc.mesozooplankton, bgc, M, O₂, T)

    GOC_fecal_pelet_production = upper_trophic_fecal_product(bgc.mesozooplankton, M, T)

    return DOC_respiration_product + DIC_respiration_product + POC_microzooplankton_mortality + GOC_mesozooplankton_mortality + GOC_fecal_pelet_production - Z_mortality - M_mortality
end
    

end # module
