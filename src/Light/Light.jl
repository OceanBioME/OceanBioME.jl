"""
Light attenuation by chlorophyll as described by [Karleskind2011](@citet) (implemented as TwoBand) and [Morel1988](@citet) (as MultiBand).
"""
module Light

export TwoBandPhotosyntheticallyActiveRadiation, 
       PrescribedPhotosyntheticallyActiveRadiation, 
       MultiBandPhotosyntheticallyActiveRadiation

using Adapt

using KernelAbstractions, Oceananigans.Units
using Oceananigans.Architectures: device, architecture, on_architecture
using Oceananigans.Utils: launch!
using Oceananigans: Center, Face, fields
using Oceananigans.Grids: node, znodes, znode
using Oceananigans.Fields: CenterField, TracerFields, location
using Oceananigans.BoundaryConditions: fill_halo_regions!, 
                                       ValueBoundaryCondition, 
                                       FieldBoundaryConditions, 
                                       regularize_field_boundary_conditions, 
                                       ContinuousBoundaryFunction,
                                       z_boundary_node,
                                       domain_boundary_indices,
                                       RightBoundary
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using OceanBioME: chlorophyll

import Adapt: adapt_structure, adapt
import Base: show, summary

import Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields, update_biogeochemical_state!, required_biogeochemical_auxiliary_fields
import Oceananigans.BoundaryConditions: _fill_top_halo!

include("2band.jl")
include("multi_band.jl")
include("prescribed.jl")

include("compute_euphotic_depth.jl")

default_surface_PAR(x, y, t) = default_surface_PAR(t)
default_surface_PAR(x_or_y, t) = default_surface_PAR(t)
default_surface_PAR(t) = 100 * max(0, cos(t * π / 12hours))

end
