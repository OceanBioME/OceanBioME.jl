module Boundaries
using Roots, Oceananigans

const defaults = (
    field_dependencies = (
        O₂ = (:O₂, ),
        CO₂ = (:DIC, :ALK)
    ),
    Sc_params = (
        O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091),
        CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219)
    ),
    β_params = (
        O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387),
        CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036)#(returns K₀ rather than β)
    ),
    pH=8.0, # initial pH value guess for air-sea flux calculation
    ρₒ = 1026, # kg m⁻³, average density at the surface of the world ocean
    conc_air = (
        CO₂=413.3,#ppmv
        O₂=209*0.0409/32#ppm from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010JC006446 to mmolO₂/m³ (https://teesing.com/en/library/tools/ppm-mg3-converter not a very good source I know)
    #this conversion is at STP(?)
    ),#may want to make these variable at some point (along with wind speed)
    molecular_weight = (
        CO₂=12+14*2#g/mol
    ),
    uₐᵥ=10#m/s https://rmets.onlinelibrary.wiley.com/doi/10.1002/joc.6957   
)

function pCO₂(DIC, ALK, T, S, params)
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    K0 = exp(-60.2409 + 9345.17/T + 23.3585*log(T/100) + S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/T - 1.468591*log(T) + (-0.138681 - 9.33291/T)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/T -1.437139*log(T) + (-0.128417 - 24.41239/T)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/T + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T) + 0.053105*sqrt(S)*T)

    ALK *= 1.e-3/params.ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC *= 1.e-3/params.ρₒ # micromoles to moles    

    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity

    H = 10^(-params.pH) # initial guess from arg list
    #potential peformance botteneck
    @inline CA(H) = ALK - (KB/(KB + H))*Boron
    @inline H_eq(H) = CA(H)*H^2 + K1*(CA(H)-DIC)*H + K1*K2*(CA(H)-2*DIC)
    H = find_zero(H_eq, H, atol=1.e-15)
    
    CO2aq = CA(H)/(K1/H + 2*K1*K2/H^2)*1e6 # Eq 11  μmol/kg
    return CO2aq/K0 #ppm
end

k(gas::Symbol, T, params)=0.39*(0.01/3600)*params.uₐᵥ^2*(Sc(T, getproperty(params.Sc_params, gas))/660)^(-0.5)#m/s, may want to add variable wind speed instead of average wind here at some point
Sc(T, params) = params.A-params.B*T+params.C*T^2-params.D*T^3

α(gas::Symbol, T, params)=β(T+273.15, S, getproperty(params.β_params, gas))/T
β(T, S, params) = exp(params.A₁+params.A₂*(100/T)+params.A₃*log(T/100)+S*(params.B₁+params.B₂*(T/100)+params.B₃*(T/100)^2))

K(gas::Symbol, T, S, params)=k(gas, T, params)*params.ρ₀*β(T+273.15, S, getproperty(params.β_params, gas)) #L=ρ\_wK₀ ->  https://reader.elsevier.com/reader/sd/pii/0304420374900152 and here K₀=β

function airseaflux(x, y, t, T::AbstractFloat, S::AbstractFloat, conc::AbstractFloat, params)
    #https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/92JC00188 could add chemical enhancement factor
    if params.gas in (:O₂, )#doing like this for flexability in the future
        return k(params.gas, T, params)*(conc-α(T, S)*getproperty(params.conc_air, params.gas))
    elseif params.gas in (:CO₂, )
        return K(params.gas, T, S, params)*(conc-getproperty(params.conc_air, params.gas))*params.ρ₀*0.001/(getproperty(params.molecular_weight, gas))#ppmv to mmolC/m³ not sure this is correct
    end
end

function airseaflux(x, y, t, T::AbstractFloat, S::AbstractFloat, DIC::AbstractFloat, ALK::AbstractFloat, params)
    if !(params.gas in (:CO₂, ))
        throw(ArgumentError("Too many arguments for gas type $(params.gas)"))
    end
    conc = pCO₂(DIC, ALK, T, S, params)
    return airseaflux(x, y, t, T, S, conc, params)
end

airseaflux(x, y, t, conc, params) = airseaflux(x, y, t, params.T(x, y, 0, t), params.S(x, y, 0, t), conc, params)
airseaflux(x, y, t, DIC, ALK, params) = airseaflux(x, y, t, params.T(x, y, 0, t), params.S(x, y, 0, t), DIC, ALK, params)

function setup(gas::Symbol; forcings=(T=nothing, S=nothing), parameters=defaults)
    #don't like this but oh well
    if !(gas in keys(parameters.Sc_params))
        throw(ArgumentError("Boundary conditions not defined for gas $gas, available are $(keys(parameters.Sc_params))"))
    end

    field_dependencies = getproperty(parameters.field_dependencies, gas)

    if (isnothing(forcings.T) & isnothing(forcings.S))
        field_dependencies = (field_dependencies..., :T, :S)
        parameters = parameters
    elseif isnothing(forcings.T)
        field_dependencies = (field_dependencies..., :T)
        parameters = merge(parameters, (S=forcings.S, ))
    elseif isnothing(forcings.S)
        field_dependencies = (field_dependencies..., :S)
        bcs_parameters = merge(parameters, (T=forcings.T, ))
    else
        field_dependencies = field_dependencies
        parameters = merge(parameters, (T=forcings.T, S=forcings.S))
    end

    parameters = merge(parameters, (gas=gas, ))

    return FluxBoundaryCondition(airseaflux, field_dependencies = field_dependencies, parameters = parameters)
end
end