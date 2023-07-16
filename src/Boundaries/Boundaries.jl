"""
Boundary conditions for air/sea and sediment flux. 

Currently implemented:

- gasexchange [Wanninkhof1992](@cite)
  - Generic air-sea flux model described by [Wanninkhof1992](@citet) but only setup for CO₂ and O₂.
  - Forces the DIC and oxygen fields, and requires temp (in centigrade) and salinity, plus
    current DIC and ALK concentration.

- Sediments
  - Soetaert [Soetaert2000](@cite)
    - simple (integrated) sediment model described by [Soetaert2000](@citet), where organic matter
       that sinks to the bottom is stored, decays into NO₃ and NH₄, and takes up O₂ in the process.
    - Extended to attribute the corresponding release of DIC.
    - Forced by O₂, NO₃, NH₄ and particle concentration in bottom cell.
  - Instant remineralisation 
    - simple model from [Aumont2015](@citet), where sinking organic matter is instantly remineralised 
      and returned to the bottom cell
    - some fraction is stored permanently in the sediment at an efficiency given by [RemineralisationFraction](@citet)
"""
module Boundaries

export Sediments, GasExchange

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units
using Oceananigans.Architectures: device
using Oceananigans.Operators: volume, Az
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Field, fractional_indices
using Oceananigans.Grids: znodes

import Adapt: adapt_structure, adapt
import Base: show, summary

include("gasexchange.jl")
include("Sediments/Sediments.jl")

end #module
