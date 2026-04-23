# convenience constructors

"""
    LOBSTER(grid; 

            nutrients = NitrateAmmonia(),
            biology = PhytoZoo(),
            detritus = TwoParticleAndDissolved(grid),
            carbonate_system = nothing,
            oxygen = nothing,

            surface_photosynthetically_active_radiation = default_surface_PAR,

            light_attenuation =
                  TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                             surface_PAR = surface_photosynthetically_active_radiation),
              
            sediment = nothing,

            scale_negatives = false,
            invalid_fill_value = NaN,

            particles = nothing,
            modifiers = nothing)

Construct an instance of the default [LOBSTER](@ref LOBSTER) biogeochemical model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to configure sinking speeds
- `nutrients`: nutrients component, defaults to `NitrateAmmonia` which has nitrate (`NO₃`) and ammonia (`NH₄`),
currently `NitrateAmmoniaIron` is the alternative which includes iron (`Fe`)
- `biology`: biological (living) component, defaults to `PhytoZoo` which has phytoplankton (`P`) and zooplankton (`Z`)
- `detritus`: non-living organic component, default to `TwoParticlesAndDissolved` which has small and large 
particles (`sPOM` and `bPOM`) and dissolved organic (`DOM`), alternatively can be `VariableRedfieldDetritus` which has 
`sPON`, `bPON`, `DON`, `sPOC`, `bPOC`, and `DOC`
- `carbonate_system`: an optional model component for the inorganic carbon components, defaults to `nothing` but
can be `CarbonateSystem` which has dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`)
- `oxygen`: an optional model component for oxygen, defaults to `nothing` but can be `Oxygen` which has oxygen (`O₂`)
- `surface_photosynthetically_active_radiation`: function for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `BiogeochemicalSediment`
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> model = LOBSTER(; grid)
LOBSTER model (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM) 
 Light attenuation: Two-band light attenuation model (Float64)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing

```
"""
LOBSTER(grid;
        nutrients = NitrateAmmonia(),
        biology = PhytoZoo(),
        detritus = TwoParticleAndDissolved(grid),
        carbonate_system = nothing,
        oxygen = nothing,

        surface_photosynthetically_active_radiation = default_surface_PAR,

        light_attenuation =
              TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                         surface_PAR = surface_photosynthetically_active_radiation),
          
        sediment = nothing,

        scale_negatives = false,
        invalid_fill_value = NaN,

        particles = nothing,
        modifiers = nothing) =
    BiologyNutrientDetritus(grid; 
                            nutrients,
                            biology,
                            detritus,
                            carbonate_system,
                            oxygen,
                            light_attenuation,
                            sediment,
                            scale_negatives,
                            invalid_fill_value,
                            particles,
                            modifiers)


"""
    NPZD(grid;
        nutrients = Nutrient(),
        biology = PhytoZoo(grid;
                            nitrate_half_saturation = 2.3868,
                            phytoplankton_maximum_growth_rate = 0.6989 / day,
                            phytoplankton_exudation_fraction = zero(grid),
                            temperature_coefficient = 1.88,
                            phytoplankton_mortality_formulation = Linear(),
                            phytoplankton_mortality_rate = (0.066 + 0.0101)/day,
                            preference_for_phytoplankton = one(grid),
                            grazing_concentration_formulation = Quadratic(),
                            grazing_half_saturation = 0.5573,
                            zooplankton_mortality_rate = 0.3395 / day,
                            zooplankton_excretion_rate = 0.0102 / day,
                            zooplankton_assimilation_fraction = 0.9116,
                            phytoplankton_sinking_speed = 0.2551/day,
                            excretion_inorganic_fraction = one(grid),
                            phytoplankton_solid_waste_fraction = 0.0101 / (0.066 + 0.0101),
                            maximum_grazing_rate = 2.1522 / day),
        detritus = Detritus(grid),

        carbonate_system = nothing,
        oxygen = nothing,

        surface_photosynthetically_active_radiation = default_surface_PAR,

        light_attenuation =
            TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                        surface_PAR = surface_photosynthetically_active_radiation),
        
        sediment = nothing,

        scale_negatives = false,
        invalid_fill_value = NaN,

        particles = nothing,
        modifiers = nothing)

Construct an instance of the default [LOBSTER](@ref LOBSTER) biogeochemical model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to configure sinking speeds
- `nutrients`: nutrients component, defaults to `NitrateAmmonia` which has nitrate (`NO₃`) and ammonia (`NH₄`),
currently `NitrateAmmoniaIron` is the alternative which includes iron (`Fe`)
- `biology`: biological (living) component, defaults to `PhytoZoo` which has phytoplankton (`P`) and zooplankton (`Z`)
- `detritus`: non-living organic component, default to `TwoParticlesAndDissolved` which has small and large 
particles (`sPOM` and `bPOM`) and dissolved organic (`DOM`), alternatively can be `VariableRedfieldDetritus` which has 
`sPON`, `bPON`, `DON`, `sPOC`, `bPOC`, and `DOC`
- `carbonate_system`: an optional model component for the inorganic carbon components, defaults to `nothing` but
can be `CarbonateSystem` which has dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`)
- `oxygen`: an optional model component for oxygen, defaults to `nothing` but can be `Oxygen` which has oxygen (`O₂`)
- `surface_photosynthetically_active_radiation`: function for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `BiogeochemicalSediment`
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> model = LOBSTER(; grid)
LOBSTER model (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM) 
 Light attenuation: Two-band light attenuation model (Float64)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing

```
"""
NPZD(grid;
     nutrients = Nutrient(),
     biology = PhytoZoo(grid;
                        nitrate_half_saturation = 2.3868,
                        phytoplankton_maximum_growth_rate = 0.6989 / day,
                        phytoplankton_exudation_fraction = zero(grid),
                        temperature_coefficient = 1.88,
                        phytoplankton_mortality_formulation = Linear(),
                        phytoplankton_mortality_rate = (0.066 + 0.0101)/day,
                        preference_for_phytoplankton = one(grid),
                        grazing_concentration_formulation = Quadratic(),
                        grazing_half_saturation = 0.5573,
                        zooplankton_mortality_rate = 0.3395 / day,
                        zooplankton_excretion_rate = 0.0102 / day,
                        zooplankton_assimilation_fraction = 0.9116,
                        phytoplankton_sinking_speed = 0.2551/day,
                        excretion_inorganic_fraction = one(grid),
                        phytoplankton_solid_waste_fraction = 0.0101 / (0.066 + 0.0101),
                        maximum_grazing_rate = 2.1522 / day,
                        light_limitation = AnalyticalLightLimitation(),
                        light_half_saturation = (0.6989/day)/(0.1953/day)),
     detritus = Detritus(grid),

     carbonate_system = nothing,
     oxygen = nothing,

     surface_photosynthetically_active_radiation = default_surface_PAR,

     light_attenuation =
           TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                      surface_PAR = surface_photosynthetically_active_radiation),
     
     sediment = nothing,

     scale_negatives = false,
     invalid_fill_value = NaN,

     particles = nothing,
     modifiers = nothing) =
    BiologyNutrientDetritus(grid; 
                            nutrients,
                            biology,
                            detritus,
                            carbonate_system,
                            oxygen,
                            light_attenuation,
                            sediment,
                            scale_negatives,
                            invalid_fill_value,
                            particles,
                            modifiers)
