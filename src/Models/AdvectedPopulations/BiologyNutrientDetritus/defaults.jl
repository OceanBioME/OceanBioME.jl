# convenience constructors

"""
    LOBSTER(; grid

              nutrients = NitrateAmmonia(),
              biology = LOBSTERPhytoZoo(),
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
- `nutrients`: nutrint component, defaults to `NitrateAmmonia` which has nitrate (`NO₃`) and ammonia (`NH₄`),
currently `NitrateAmmoniaIron` is the alternative which includes iron (`Fe`)
- `biology`: biological (living) component, defaults to `LOBSTERPhytoZoo` which has phytoplankton (`P`) and zooplankton (`Z`)
- `detritus`: non-living organic component, default to `TwoParticlesAndDissolved` which has small and large 
particles (`sPOM` and `bPOM`) and dissolved organic (`DOM`), alternativly can be `VariableRedfieldDetritus` which has 
`sPON`, `bPON`, `DON`, `sPOC`, `bPOC`, and `DOC`
- `carbonate_system`: an optional model component for the inorganic carbon components, defaults to `nothing` but
can be `CarbonateSystem` which has dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`)
- `oxygen`: an optional model component for oxygen, defaults to `nothing` but can be `Oxygen` which has oxygen (`O₂`)
- `surface_photosynthetically_active_radiation`: funciton for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
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
function LOBSTER(; grid, 

                   nutrients = NitrateAmmonia(),
                   biology = LOBSTERPhytoZoo(),
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

    lobster = BiologyNutrientDetritus(nutrients, biology, detritus, carbonate_system, oxygen)

    if scale_negatives
        scaler = ScaleNegativeTracers(lobster, grid; invalid_fill_value)
        if isnothing(modifiers)
            modifiers = scaler
        elseif modifiers isa Tuple
            modifiers = (modifiers..., scaler)
        else
            modifiers = (modifiers, scaler)
        end
    end
    
    return Biogeochemistry(lobster;
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end