"
Boundary conditions for air/sea and sediment flux. 

Currently implimented:
- airseaflux
    - Generic air sea flux  model describted by Wanninkhof, 1992 but only setup for CO₂ and O₂
    - Forces the DIC field, and requires temp (in centigrade) and salinity, plus current DIC and ALK concentration
    - Extended to provide boundary conditon on DIC at redfield ratio (mol DIC/mol NH₄) of 106/16 from Table 4 caption
- sediment
    - simple (integrated) sediment model described by Soetaert, Middelburg, Herman and Buis, 2000 
    where organic matter (D and DD) that sinks to the bottom is stored and decays into NO₃ and NH₄, 
    and takes up O₂ in the process. 
    - Extended to attribute the corresponding release of DIC
    - Forced by O₂, D, and DD concentration in bottom cell

References:
Soetaert, K., Middelburg, J., Herman, P. and Buis, K., 2000. On the coupling of benthic and pelagic biogeochemical models. Earth-Science Reviews, 51(1-4), pp.173-201.
Wanninkhof, R., 1992. Relationship between wind speed and gas exchange over the ocean. Journal of Geophysical Research, 97(C5), p.7373.
"
module Boundaries
using Roots, Oceananigans, KernelAbstractions, Adapt
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device

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
        uₐᵥ=10#m/s https://rmets.onlinelibrary.wiley.com/doi/10.1002/joc.6957
    ),
    sediment = (
        #https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.1996.41.8.1651
        λᵣᵣ = 2/year,# 1/year to 1/s
        λᵣ = 0.2/year,#s   
        Rdᵣᵣ = 0.1509,#mmol N/mmol C
        Rdᵣ = 0.13,#mmol N/mmol C
        Rd_red = 106/16#mmol C/mmol N
    )
)

include("airseaflux.jl")
include("sediment.jl")
end
