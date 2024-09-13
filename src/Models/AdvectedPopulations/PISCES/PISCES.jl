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

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
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
(bgc::PISCES)(val_name::carbon_system, args...) = bgc.carbon_system(val_name, bgc, args...)


@inline biogeochemical_auxiliary_fields(bgc::PISCES) = 
    (zₘₓₗ = bgc.mixed_layer_depth, 
     zₑᵤ = bgc.euphotic_depth, 
     Si̅ = bgc.yearly_maximum_silicate, 
     D_dust = bgc.dust_deposition, 
     Ω = bgc.calcite_saturation,
     κ = bgc.mean_mixed_layer_vertical_diffusivity,
     mixed_layer_PAR = bgc.mean_mixed_layer_light)

@inline required_biogeochemical_tracers(::PISCES) = 
    (:P, :D, :Z, :M, :Pᶜʰˡ, :Dᶜʰˡ, :Pᶠᵉ, :Dᶠᵉ, :Dˢⁱ, 
     :DOC, :POC, :GOC, :SFe, :BFe, :PSi, 
     :NO₃, :NH₄, :PO₄, :Fe, :Si, 
     :CaCO₃, :DIC, :Alk, :O₂, :T)

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
include("update_state.jl")

function PISCES(; grid,
                  nanophytoplankton = 
                    Phytoplankton(growth_rate = GrowthRespirationLimitedProduction(dark_tollerance = 3days),
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
                    Phytoplankton(growth_rate = GrowthRespirationLimitedProduction(dark_tollerance = 4days),
                                  nutrient_limitation = 
                                    NitrogenIronPhosphateSilicateLimitation(minimum_ammonium_half_saturation = 0.039,
                                                                            minimum_nitrate_half_saturation = 0.39, 
                                                                            minimum_phosphate_half_saturation = 2.4,
                                                                            half_saturation_for_iron_uptake = 3.0,
                                                                            silicate_limited = true),
                                  blue_light_absorption = 1.6, 
                                  green_light_absorption = 0.69, 
                                  red_light_absorption = 0.7,
                                  maximum_quadratic_mortality = 0.03,
                                  maximum_chlorophyll_ratio = 0.05),
                  
                  microzooplankton = Zooplankton(maximum_grazing_rate = 3/day,
                                                 preference_for_nanophytoplankton = 1.0,
                                                 preference_for_diatoms = 0.5,
                                                 preference_for_particulates = 0.1,
                                                 preference_for_microzooplankton = 0.0
                                                 quadratic_mortality = 0.004,
                                                 linear_mortality = 0.03,
                                                 maximum_growth_efficiency = 0.3,
                                                 maximum_flux_feeding_rate = 0.0,
                                                 undissolved_calcite_fraction = 0.5),

                  mesozooplankton = Zooplankton(maximum_grazing_rate = 0.75/day,
                                                preference_for_nanophytoplankton = 0.3,
                                                preference_for_diatoms = 1.0,
                                                preference_for_particulates = 0.3,
                                                preference_for_microzooplankton = 1.0,
                                                quadratic_mortality = 0.03,
                                                linear_mortality = 0.005,
                                                maximum_growth_efficiency = 0.35,
                                                maximum_flux_feeding_rate = 2.0e-3,
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
                  
                  mixed_layer_shear = 1.0,
                  background_shear = 0.01, 
                  
                  latitude = PrescribedLatitude(45),
                  day_length,
                  
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
                                        nitrogen_redfield_ratio, phosphate_redfield_ratio,
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

end # module
