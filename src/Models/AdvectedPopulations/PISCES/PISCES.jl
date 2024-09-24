"""
Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model.

This is *not* currently an official version supported by the PISCES community
and is not yet verified to be capable of producing results mathcing that of the 
operational PISCES configuration. This is a work in progress, please open an 
issue or discusison if you'd like to know more.

Notes to developers
===================
Part of the vision for this implementation of PISCES is to harness the features
of Julia that would allow it to be fully modular. An obvious step to improve the
ease of this would be to do some minor refactoring to group the phytoplankton 
classes, and zooplankton classes together, and for the other groups to generically 
call the whole lot. This may cause some issues with argument passing, and although
it may not be the best way todo it my first thought is to pass them round as named
tuples built from something like,
```
phytoplankton_tracers = phytoplankton_arguments(bgc.phytoplankton, args...)
```

"""
module PISCESModel

export PISCES, DepthDependantSinkingSpeed, PrescribedLatitude, ModelLatitude

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

import Base: show, summary

struct PISCES{PP, ZP, DM, PM, NI, FE, SI, OX, PO, CA, CE, FT, LA, DL, ML, EU, MS, VD, MP, CC, CS, SS} <: AbstractBiogeochemistry
                            phytoplankton :: PP

                              zooplankton :: ZP

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
                     silicate_climatology :: MS

    mean_mixed_layer_vertical_diffusivity :: VD
                   mean_mixed_layer_light :: MP

                         carbon_chemistry :: CC
                       calcite_saturation :: CS

                       sinking_velocities :: SS
end

@inline required_biogeochemical_tracers(bgc::PISCES) = 
    (required_biogeochemical_tracers(bgc.zooplankton)...,
     required_biogeochemical_tracers(bgc.phytoplankon)...,
     required_biogeochemical_tracers(bgc.dissolved_organic_matter)...,
     required_biogeochemical_tracers(bgc.particulate_organic_matter)...,
     required_biogeochemical_tracers(bgc.nitrogen)...,
     required_biogeochemical_tracers(bgc.phosphate)...,
     required_biogeochemical_tracers(bgc.iron)...,
     required_biogeochemical_tracers(bgc.silicate),
     required_biogeochemical_tracers(bgc.carbon_system)...
     required_biogeochemical_tracers(bgc.oxygen)...,
     :T, :S)

@inline required_biogeochemical_auxiliary_fields(::PISCES) =
    (:zₘₓₗ, :zₑᵤ, :Si′, :Ω, :κ, :mixed_layer_PAR, :wPOC, :wGOC, :PAR, :PAR₁, :PAR₂, :PAR₃)

@inline biogeochemical_auxiliary_fields(bgc::PISCES) = 
    (zₘₓₗ = bgc.mixed_layer_depth, 
     zₑᵤ = bgc.euphotic_depth, 
     Si′ = bgc.silicate_climatology, 
     Ω = bgc.calcite_saturation,
     κ = bgc.mean_mixed_layer_vertical_diffusivity,
     mixed_layer_PAR = bgc.mean_mixed_layer_light,
     wPOC = bgc.sinking_velocities.POC,
     wGOC = bgc.sinking_velocities.GOC)

biogeochemical_drift_velocity(bgc::PISCES, val_name) = 
    biogeochemical_drift_velocity(bgc.particulate_organic_matter, val_name)

include("zooplankton/zooplankton.jl")

using .Zooplankton

include("phytoplankton/phytoplankton.jl")

using .Phytoplankton

include("dissolved_organic_matter/dissolved_organic_matter.jl")

using .DissolvedOrganicMatter

include("particulate_organic_matter/particulate_organic_matter.jl")

using .ParticulateOrganicMatter

include("nitrogen/nitrogen.jl")

using .Nitrogen

include("iron/iron.jl")

using .Iron

include("silicate.jl")

using .Silicates

include("oxygen.jl")

using .OxygenModels

include("phosphate.jl")

using .Phosphates

include("inorganic_carbon.jl")

using .InorganicCarbons

"""
    PISCES(; grid,
             

             # from Aumount 2005 rather than 2015 since it doesn't work the other way around
             first_anoxia_thresehold = 6.0,
             second_anoxia_thresehold = 1.0,

             nitrogen_redfield_ratio = 16/122,
             phosphate_redfield_ratio = 1/122,
             
             mixed_layer_shear = 1.0,
             background_shear = 0.01, 
             
             latitude = PrescribedLatitude(45),
             day_length = day_length_function,
             
             mixed_layer_depth = Field{Center, Center, Nothing}(grid),
             euphotic_depth = Field{Center, Center, Nothing}(grid),

             silicate_climatology = ConstantField(7.5),

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
                               GOC = Field(KernelFunctionOperation{Center, Center, Face}(DepthDependantSinkingSpeed(), 
                                                                                         grid, 
                                                                                         mixed_layer_depth, 
                                                                                         euphotic_depth)),
             open_bottom = true,

             scale_negatives = false,
             invalid_fill_value = NaN,
             
             sediment = nothing,
             particles = nothing,
             modifiers = nothing)

Constructs an instance of the PISCES biogeochemical model.


Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `nanophytoplankton`: nanophytoplankton (`P`, `PChl`, `PFe``) evolution parameterisation such as `MixedMondoPhytoplankton`
- `diatoms`: diatom (`D`, `DChl`, `DFe`, `DSi`) evolution parameterisation such as `MixedMondoPhytoplankton`
- `microzooplankton`: microzooplankton (`Z`) evolution parameterisation
- `mesozooplankton`: mesozooplankton (`M`) evolution parameterisation
- `dissolved_organic_matter`: parameterisaion for the evolution of dissolved organic matter (`DOC`)
- `particulate_organic_matter`: parameterisation for the evolution of particulate organic matter (`POC`, `GOC`, `SFe`, `BFe`, `PSi`)
- `nitrogen`: parameterisation for the nitrogen compartements (`NH₄` and `NO₃`)
- `iron`: parameterisation for iron (`Fe`), currently the "complex chemistry" of Aumount 2015 is not implemented
- `silicate`: parameterisaion for silicate (`Si`)
- `oxygen`: parameterisaion for oxygen (`O₂`)
- `phosphate`: parameterisaion for phosphate (`PO₄`)
- `calcite`: parameterisaion for calcite (`CaCO₃`)
- `carbon_system`: parameterisation for the evolution of the carbon system (`DIC` and `Alk`alinity)
- `first_anoxia_thresehold` and `second_anoxia_thresehold`: thresholds in anoxia parameterisation
- `nitrogen_redfield_ratio` and `phosphate_redfield_ratio`: the assumed element ratios N/C and P/C 
- `mixed_layer_shear` and `background_shear`: the mixed layer and background shear rates, TODO: move this to a computed field
- `latitude`: model latitude, should be `PrescribedLatitude` for `RectilinearGrid`s and `ModelLatitude` for grids providing their own latitude
- `day_length`: parameterisation for day length based on time of year and latitude, you may wish to change this to (φ, t) -> 1day if you
   want to ignore the effect of day length, or something else if you're modelling a differen planet
- `mixed_layer_depth`: an `AbstractField` containing the mixed layer depth (to be computed during update state)
- `euphotic`: an `AbstractField` containing the euphotic depth, the depth where light reduces to 1/1000 of 
   the surface value (computed during update state)
- `silicate_climatology`: an `AbstractField` containing the silicate climatology which effects the diatoms silicate
   half saturation constant
- `mean_mixed_layer_vertical_diffusivity`: an `AbstractField` containing the mean mixed layer vertical diffusivity 
   (to be computed during update state)
- `mean_mixed_layer_light`: an `AbstractField` containing the mean mixed layer light (computed during update state)
- `carbon_chemistry`: the `CarbonChemistry` model used to compute the calicte saturation
- `calcite_saturation`: an `AbstractField` containing the calcite saturation  (computed during update state)
- `surface_photosynthetically_active_radiation`: funciton for the photosynthetically available radiation at the surface
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sinking_speed`: named tuple of constant sinking speeds, or fields (i.e. `ZFaceField(...)`) for any tracers which sink 
  (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
- `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

All parameterisations default to the operaitonal version of PISCES as close as possible.

Notes
=====
Currently only `MixedMondoPhytoplankton` are implemented, and some work should be done to generalise
the classes to a single `phytoplankton` if more classes are required (see 
`OceanBioME.Models.PISCESModel` docstring). Similarly, if a more generic `particulate_organic_matter`
was desired a way to specify arbitary tracers for arguments would be required.
"""
function PISCES(; grid,
                  phytoplankton = MixedMondoNanoAndDiatoms(),
                  zooplankton = MicroAndMezoZooplankton(),
                  dissolved_organic_matter = DissolvedOrganicCarbon(),
                  particulate_organic_matter = TwoCompartementCarbonIronParticles(),
                  
                  nitrogen = NitrateAmmonia(),
                  iron = SimpleIron(),
                  silicate = Silicate(),
                  oxygen = Oxygen(),
                  phosphate = Phosphate(),
                  
                  carbon_system = InorganicCarbon(),

                  # from Aumount 2005 rather than 2015 since it doesn't work the other way around
                  first_anoxia_thresehold = 6.0,
                  second_anoxia_thresehold = 1.0,

                  nitrogen_redfield_ratio = 16/122,
                  phosphate_redfield_ratio = 1/122,
                  
                  mixed_layer_shear = 1.0,
                  background_shear = 0.01, 
                  
                  latitude = PrescribedLatitude(45),
                  day_length = day_length_function,
                  
                  mixed_layer_depth = Field{Center, Center, Nothing}(grid),
                  euphotic_depth = Field{Center, Center, Nothing}(grid),

                  silicate_climatology = ConstantField(7.5),

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
                                    GOC = Field(KernelFunctionOperation{Center, Center, Face}(DepthDependantSinkingSpeed(), 
                                                                                              grid, 
                                                                                              mixed_layer_depth, 
                                                                                              euphotic_depth))),
                  open_bottom = true,

                  scale_negatives = false,
                  invalid_fill_value = NaN,
                  
                  sediment = nothing,
                  particles = nothing,
                  modifiers = nothing)

    @warn "This implementation of PISCES is in early development and has not yet been validated"

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
                                        silicate_climatology,
                                        mean_mixed_layer_vertical_diffusivity, 
                                        mean_mixed_layer_light,
                                        carbon_chemistry, calcite_saturation,
                                        sinking_velocities)

    if scale_negatives
        scalers = ScaleNegativeTracers(underlying_biogeochemistry, grid; invalid_fill_value)
        if isnothing(modifiers)
            modifiers = scalers
        elseif modifiers isa Tuple
            modifiers = (modifiers..., scalers...)
        else
            modifiers = (modifiers, scalers...)
        end
    end

    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end

end # module
