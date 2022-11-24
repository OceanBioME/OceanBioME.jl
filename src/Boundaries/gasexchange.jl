#####
##### Gas exchange model of [Wanninkhof1992](@cite)
#####
# TODO: Impliment Ho et al. 2006 wind speed dependence

##### 
##### Carbonate chemistry to determine pCO₂
#####

@inline CA_eq(H, params) = params.ALK - (params.KB/(params.KB + H))*params.Boron
@inline H_eq(H, params) = CA_eq(H, params)*H^2 + params.K1*(CA_eq(H, params)-params.DIC)*H + params.K1*params.K2*(CA_eq(H, params)-2*params.DIC)

function pCO₂(DIC, ALK, T, S, ρₒ, pH)
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    ALK *= 1.e-3/ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC *= 1.e-3/ρₒ # micromoles to moles    

    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity

    K0 = exp(-60.2409 + 9345.17/T + 23.3585*log(T/100) + S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/T - 1.468591*log(T) + (-0.138681 - 9.33291/T)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/T -1.437139*log(T) + (-0.128417 - 24.41239/T)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/T + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T) + 0.053105*sqrt(S)*T)

    H = 10^(-pH) # initial guess from arg list

    p = (DIC=DIC, ALK=ALK, K0=K0, K1=K1, K2=K2, KB=KB, Boron=Boron)

    H = find_zero(H_eq, H, atol=1e-100, p=p)
    CA = CA_eq(H, p)
    
    #pH = -log10(H)
    CO2aq = 1e6*CA/(K1/H + 2*K1*K2/H^2) # Eq 11  μmol/kg
    return CO2aq/K0 #μatm
end

k(T, uₐᵥ, Sc_params) = 0.39 * (0.01 / 3600) * uₐᵥ ^ 2 * (Sc(T, Sc_params) / 660) ^ (-0.5)# m/s, may want to add variable wind speed instead of average wind here at some point
Sc(T, params) = params.A - params.B * T + params.C * T ^ 2 - params.D * T ^ 3

α(T, S, β_params) = β(T + 273.15, S, β_params) * (T + 273.15) * 0.00367#/(T+273.15) - disagree with origional paper but this matches dimensionless Henry coefficiet of 3.2x10⁻² at 298.15K, S=0. See https://www.wikiwand.com/en/Henry%27s_law 
β(T, S, params) = exp(params.A₁ + params.A₂ * (100 / T) + params.A₃ * log(T / 100) + S * (params.B₁ + params.B₂ * (T / 100) + params.B₃ * (T / 100) ^ 2))

#now fairly sure that this is giving the correct result for CO₂  as βρ is the henrys coefficient which sould be ∼34mol/m³ atm
K(T, S, uₐᵥ, Sc_params, β_params, ρₒ) = k(T, uₐᵥ, Sc_params) * β(T+273.15, S, β_params) * ρₒ #L=ρ\_wK₀ ->  https://reader.elsevier.com/reader/sd/pii/0304420374900152 and here K₀=β

#####
##### Boundary condition setup
#####

struct GasExchange{G, ScP, βP, FT, T, S}
    gas :: G

    schmidt_params :: ScP
    solubility_params :: βP
    pH_initial_guess :: FT
    ocean_density :: FT
    air_concentration :: FT
    air_pressure :: FT
    average_wind_speed :: FT

    temperature :: T
    salinity :: S
end

"""
    GasExchange(;gas,
                schmidt_params::ScP = (CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219),
                                O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091))[gas],
                solubility_params::βP = (CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036),
                                    O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387))[gas],
                pH_initial_guess::FT = 8.0,
                ocean_density::FT = 1026, # kg/m³
                air_concentration::FT = (CO₂ = 413.4, O₂ = 9352.7)[gas], # ppmv, mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
                air_pressure::FT = 1.0, # atm
                average_wind_speed::FT = 10, # m/s
                field_dependencies = (CO₂ = (:DIC, :ALK), O₂ = (:OXY, ))[gas],
                temperature::T = nothing,
                salinity::S = nothing)

Constructs an Oceananigans `FluxBoundaryCondition` for the exchange of `gas` with the relivant tracer (i.e. DIC for CO₂ and oxygen for O₂).
Please see note for other gases.

Keyword arguments
=================

    - `gas`: (required) the gas to be exchanged, if `:CO₂` or `:O₂` are specified then all other settings may be infered
    - `schmidt_params` : named tuple of parameters for calculating the Schmidt number using the parameterisation of [Wanninkhof1992](@cite)
    - `solubility_params` : named tuple of parameters for calculating the solubility (for O₂ the Bunsen solubility and CO₂ K₀, see note)
    - `pH_initial_guess` : initial guess of pH for calculating pCO₂ - is not used for other gases
    - `ocean_density` : density of the ocean in kg/m³
    - `air_concentratio` : concentration of the gas in air in relivant units (i.e. ppmv for CO₂ and mmol O₂/m³ for O₂)
    - `air_pressure` : air pressure in atm (only used for CO₂)
    - `average_wind_speed` : average wind speed at 10m used to calculate the gas transfer velocity by the [Wanninkhof1992](@cite) parameterisation
    - `field_dependencies` : tracer fields that gas exchange depends on, if the defaults have different names in your model you can specify as long as they are in the same order
    - `temperature` : either `nothing` to track a temperature tracer field, or a function or shape `f(x, y, z, t)` for the temperature in °C
    - `salinity` : either `nothing` to track a salinity tracer field, or a function or shape `f(x, y, z, t)` for the salinity in ‰

Note
=====

This model is fully capable of exchanging any gas but the parameters have only been configured for CO₂ and O₂, and the specific formulaiton
is only ensured for these gasses. For any gas where the [Wanninkhof1992](@cite) parameterisation returns the Bunsen Solubility Coefficient
this model will work out of the box and can just be passed new parameters. For the other solubility types (i.e. K₀, K' and f) you will need
to overload the `(gasexchange::GasExchange)` function to ensure the correct formulaiton.
"""

function GasExchange(;gas,
                      schmidt_params::ScP = (CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219),
                                      O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091))[gas],
                      solubility_params::βP = (CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036),
                                            O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387))[gas],
                      pH_initial_guess::FT = 8.0,
                      ocean_density::FT = 1026.0, # kg/m³
                      air_concentration::FT = (CO₂ = 413.4, O₂ = 9352.7)[gas], # ppmv, mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
                      air_pressure::FT = 1.0, # atm
                      average_wind_speed::FT = 10.0, # m/s
                      field_dependencies = (CO₂ = (:DIC, :ALK), O₂ = (:OXY, ))[gas],
                      temperature::T = nothing,
                      salinity::S = nothing) where {ScP, βP, FT, T, S}

    gas = Val(gas)
    G = typeof(gas)

    gasexchange =  GasExchange{G, ScP, βP, FT, T, S}(gas, 
                                                     schmidt_params, 
                                                     solubility_params, 
                                                     pH_initial_guess, 
                                                     ocean_density, 
                                                     air_concentration, 
                                                     air_pressure, 
                                                     average_wind_speed, 
                                                     temperature, 
                                                     salinity)

    if isnothing(temperature)
        field_dependencies = (field_dependencies..., :T)
    end

    if isnothing(salinity)
        field_dependencies = (field_dependencies..., :S)
    end

    return FluxBoundaryCondition(gasexchange_function, field_dependencies = field_dependencies, parameters = gasexchange)
end

# Hack this nicer format into Oceananians boundary conditions because `typeof(gasexhcange) = DataType` and `BoundaryCondition` has a special use for `DataType` in the first argument so doesn't work
@inline gasexchange_function(x, y, t, args...) = @inbounds args[end](x, y, t, args[1:end-1]...)

@inline (gasexchange::GasExchange)(x, y, t, conc) = gasexchange(x, y, t, conc, gasexchange.temperature(x, y, 0.0, t), gasexchange.salinity(x, y, 0.0, t))
@inline (gasexchange::GasExchange)(x, y, t, DIC, ALK) = gasexchange(x, y, t, DIC, ALK, gasexchange.temperature(x, y, 0.0, t), gasexchange.salinity(x, y, 0.0, t))
@inline function (gasexchange::GasExchange)(x, y, t, DIC, ALK, T, S) 
    conc = pCO₂(DIC, ALK, T + 273.15, S, gasexchange.ocean_density, gasexchange.pH_initial_guess)
    return gasexchange(x, y, t, conc, T, S)
end
@inline function (gasexchange::GasExchange)(x, y, t, conc, T, S) 
    return k(T, gasexchange.average_wind_speed, gasexchange.schmidt_params) * (conc - α(T, S, gasexchange.solubility_params) * gasexchange.air_concentration)
end
@inline function (gasexchange::GasExchange{<:Val{:CO₂}, <:Any, <:Any, <:Any, <:Any, <:Any})(x, y, t, conc, T, S) 
    return K(T, S, gasexchange.average_wind_speed, gasexchange.schmidt_params, gasexchange.solubility_params, gasexchange.ocean_density) * (conc - gasexchange.air_concentration * gasexchange.air_pressure) / 1000#μmol/m²s to mmolC/m²s not sure this is correct
end