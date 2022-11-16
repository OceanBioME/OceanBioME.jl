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

export Sediments

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device
using Oceananigans.Operators: volume, Az

const defaults = (
    airseaflux = (
        field_dependencies = (
            O₂ = (:OXY, ),
            CO₂ = (:DIC, :ALK)
        ),
        Sc_params = (
            O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091),
            CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219)
        ),
        β_params = (
            O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387),
            CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036)#(returns K₀ rathe1 than β)
        ),
        pH=8.0, # initial pH value guess for air-sea flux calculation
        ρₒ = 1026, # kg m⁻³, average density at the surface of the world ocean
        conc_air = (
            CO₂= 413.3,#ppmv
            O₂= 9352.7#mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
        #this conversion is at STP(?)
        ),#may want to make these variable at some point (along with wind speed)
        pAir = 1.0,# atm
        uₐᵥ=10#m/s https://rmets.onlinelibrary.wiley.com/doi/10.1002/joc.6957
    ),
)

include("gasexchange.jl")
include("Sediments/Sediments.jl")
end