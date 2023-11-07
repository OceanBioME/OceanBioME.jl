"""
Light attenuation by chlorophyll as described by [Karleskind2011](@citet) (implemented as twoBand) and [Morel1988](@citet).
"""
module Light

export TwoBandPhotosyntheticallyActiveRadiation, update_PAR!

using KernelAbstractions, Oceananigans.Units
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans: Center, Face, fields
using Oceananigans.Grids: node, znodes
using Oceananigans.Fields: CenterField, TracerFields
using Oceananigans.BoundaryConditions: fill_halo_regions!, 
                                       ValueBoundaryCondition, 
                                       FieldBoundaryConditions, 
                                       regularize_field_boundary_conditions, 
                                       ContinuousBoundaryFunction,
                                       z_boundary_node,
                                       domain_boundary_indices,
                                       RightBoundary
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

import Adapt: adapt_structure, adapt
import Base: show, summary

import Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields, update_biogeochemical_state!, required_biogeochemical_auxiliary_fields
import Oceananigans.BoundaryConditions: _fill_top_halo!

include("2band.jl")
include("morel.jl")

end
