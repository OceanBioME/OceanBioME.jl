"""
Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model.

Tracers
=======
# see the others for the formatting here...you might also want to change some of the units buit molC/L is the origional (described at the start of sec 4)
* Nano-phytoplankton: P (mol C/L)

Required submodels
==================
# you will need something like this, they use a different PAR model but I wouldn't worry about that for now, you might also need temperatre and salinity (not sure)
* Photosynthetically available radiation: PAR (W/m²)

"""
module PISCESModel

export PISCES

using Oceananigans.Units
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR
using OceanBioME: setup_velocity_fields, show_sinking_velocities, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.BoxModels: BoxModel
using OceanBioME.Boundaries.Sediments: sinking_flux

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

import OceanBioME: redfield, conserved_tracers

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt
import Base: show, summary

import OceanBioME.Boundaries.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

struct PISCES{FT, W} <: AbstractContinuousFormBiogeochemistry
    parameter_1 :: FT # add list of parameters here, assuming theyre all just numbers FT will be fine for advect_particles_kernel

    sinking_velocities :: W

    function PISCES(parameter_1::FT, # then do the same here (this is all just annoying boiler plate but we need it to make the next function work)
                    sinking_velocities::W) where {FT, W}

        return new{FT, W}(parameter_1, # and list them all again here...
                          sinking_velocities)
    end
end

"""
    PISCES(; grid,
             parameter_1::FT = 1.0, # now you can finally put the values here

             surface_photosynthetically_active_radiation = default_surface_PAR,

             light_attenuation_model::LA =
                 TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                                     surface_PAR = surface_photosynthetically_active_radiation),

             # just keep all this stuff for now but you can ignore it
             sediment_model::S = nothing,

             sinking_speeds = (sPOM = 3.47e-5, bPOM = 200/day),
             open_bottom::Bool = true,

             scale_negatives = false,

             particles::P = nothing,
             modifiers::M = nothing)

Construct an instance of the [PISCES](@ref PISCES) biogeochemical model. 

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to calculate sinking
- `parameter_1`...: PISCES parameter values
- `surface_photosynthetically_active_radiation`: funciton (or array in the future) for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `sinking_speed`: named tuple of constant sinking, of fields (i.e. `ZFaceField(...)`) for any tracers which sink (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
- `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME

julia> using Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> model = PISCES(; grid)
PISCES{Float64} ... # we can fix this later
```
"""
function PISCES(; grid, # finally the function
                  parameter_1::FT = 1.0, # now you can finally put the values here

                  surface_photosynthetically_active_radiation = default_surface_PAR,

                  light_attenuation_model::LA =
                    TwoBandPhotosyntheticallyActiveRadiation(; grid, 
                                            surface_PAR = surface_photosynthetically_active_radiation),

                  # just keep all this stuff for now but you can ignore it
                  sediment_model::S = nothing,

                  sinking_speeds = (sPOM = 3.47e-5, bPOM = 200/day),
                  open_bottom::Bool = true,

                  scale_negatives = false,

                  particles::P = nothing,
                  modifiers::M = nothing) where {FT, LA, S, P, M}

    if !isnothing(sediment_model) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    underlying_biogeochemistry = PISCES(parameter_1, # last time to put them all in

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
                           light_attenuation = light_attenuation_model, 
                           sediment = sediment_model, 
                           particles,
                           modifiers)
end

@inline required_biogeochemical_tracers(::PISCES) = (:P, ) # list all the parameters here, also if you need T and S put them here too

@inline required_biogeochemical_auxiliary_fields(::PISCES) = (:PAR, )

# for sinking things like POM this is how we tell oceananigans ther sinking speed
@inline function biogeochemical_drift_velocity(bgc::PISCES, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end

# don't worry about this for now
adapt_structure(to, pisces::PISCES) =
    PISCES(adapt(to, pisces.parameter_1),
           adapt(to, pisces.sinking_velocities))

# you can updatye these if you want it to have a pretty way of showing uyou its a pisces model
summary(::PISCES{FT}) where {FT} = string("PISCES{$FT}") 

show(io::IO, model::PISCES) where {FT, B, W}  = print(io, string("Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model")) # maybe add some more info here

@inline maximum_sinking_velocity(bgc::PISCES) = maximum(abs, bgc.sinking_velocities.bPOM.w) # might need ot update this for wghatever the fastest sinking pareticles are

# write most of the code here (i.e. make a file falled phytoplankton.jl and then include it here)
include("phytoplankton.jl")

# to work with the sediment model we need to tell in the redfield ratio etc. of some things, but for now we can ignore
@inline redfield(i, j, k, val_tracer_name, bgc::PISCES, tracers) = NaN

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline remineralisation_receiver(::PISCES) = :NH₄

# this is for positivity preservation, if you can work it out it would be great, I don't think PISCES conserves C but probably does Nitrogen
@inline conserved_tracers(::PISCES) = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

@inline sinking_tracers(::PISCES) = (:sPOM, :bPOM) # please list them here
end # module
