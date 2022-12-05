"
Boundary conditions for air/sea and sediment flux. 

Currently implimented:
- gasexchange ([Wanninkhof1992](@cite))
    - Generic air sea flux  model describted by Wanninkhof, 1992 but only setup for CO₂ and O₂
    - Forces the DIC and ocygen fields, and requires temp (in centigrade) and salinity, plus current DIC and ALK concentration
- Sediments
    - Soetaert ([Soetaert2000](@cite))
        - simple (integrated) sediment model described by Soetaert, Middelburg, Herman and Buis, 2000 
        where organic matter (D and DD) that sinks to the bottom is stored and decays into NO₃ and NH₄, 
        and takes up O₂ in the process. 
        - Extended to attribute the corresponding release of DIC
        - Forced by O₂, NO₃, NH₄ and particle concentration in bottom cell
"
module Boundaries

export Sediments, GasExchange

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device
using Oceananigans.Operators: volume, Az
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Field, interpolate

include("gasexchange.jl")
include("Sediments/Sediments.jl")
end