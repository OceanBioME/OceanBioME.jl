"""
MARBLModel

Tracers:
- Z
- N
- DIC
- Alk
- O₂
- NO₃
- NH₄
- Fe
- DOM
"""
module MARBLModel

export MARBL

using Oceananigans.Units

using Oceananigans: KernelFunctionOperation
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField, ConstantField, Center, Face

using OceanBioME.Light: MultiBandPhotosyntheticallyActiveRadiation, default_surface_PAR, compute_euphotic_depth!
using OceanBioME: setup_velocity_fields, show_sinking_velocities, Biogeochemistry, DiscreteBiogeochemistry, ScaleNegativeTracers, CBMDayLength
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import OceanBioME: redfield, conserved_tracers, maximum_sinking_velocity, chlorophyll

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!

import Base: show, summary

struct MARBL{PP, ZP, NC, DO, PO, CC, CS, SS} <: AbstractBiogeochemistry
       phytoplankton :: PP
         zooplankton :: ZP
            nitrogen :: NC # not sure how best to modularise but I can envisage making it so we can fix the stochiometry and remove some elements so splitting elementally might be best

  dissolved_organics :: DO
particulate_organics :: PO

    carbon_chemistry :: CC
  calcite_saturation :: CS # I think this is needed, see PISCES

  sinking_velocities :: SS
end

@inline required_biogeochemical_tracers(bgc::MARBL) = 
    (required_biogeochemical_tracers(bgc.phytoplankton)...,
     required_biogeochemical_tracers(bgc.zooplankton)...,
     required_biogeochemical_tracers(bgc.dissolved_organics)...,
     required_biogeochemical_tracers(bgc.particulate_organics)...,
     required_biogeochemical_tracers(bgc.nitrogen)...,
     :T, :S)

@inline required_biogeochemical_auxiliary_fields(::MARBL) = (:Ω, :PAR) # only one PAR band I think?

@inline biogeochemical_auxiliary_fields(bgc::MARBL) = (; Ω = bgc.calcite_saturation)

(bgc::MARBL)(i, j, k, grid, val_name, clock, fields, auxiliary_fields) = zero(grid)

(bgc::DiscreteBiogeochemistry{<:MARBL})(i, j, k, grid, val_name, clock, fields) =
    bgc.underlying_biogeochemistry(i, j, k, grid, val_name, clock, fields, biogeochemical_auxiliary_fields(bgc))

include("common.jl")
include("update_state.jl")

include("phytoplankton/phytoplankton.jl")

using .Phytoplankton

# add the other groups here

"""
    MARBL(; grid,
            ...)

Constructor
"""
function MARBL(; grid,
                  
                  # whatever else needed here...

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

    @warn "This implementation of MARBL is in early development and has not yet been validated against the operational version"

    if !isnothing(sediment) && !open_bottom
        @warn "You have specified a sediment model but not `open_bottom` which will not work as the tracer will settle in the bottom cell"
    end

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    sinking_velocities = merge(sinking_velocities, (; grid)) # we need to interpolate the fields so we need this for flux feeding inside a kernel - this might cause problems...

    underlying_biogeochemistry = MARBL(phytoplankton, ) # the other things in here

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

    # put it all in a Biogeochemistry
    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation, 
                           sediment, 
                           particles,
                           modifiers)
end

end # module
