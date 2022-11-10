#####
##### Gas exchange model of [Wanninkhof1992](@cite)
#####
# TODO: Impliment Ho et al. 2006 wind speed dependence

##### 
##### Carbonate chemistry to determine pCO₂
#####

@inline CA_eq(H, params) = params.ALK - (params.KB/(params.KB + H))*params.Boron

@inline H_eq(H, params) = CA_eq(H, params)*H^2 + params.K1*(CA_eq(H, params)-params.DIC)*H + params.K1*params.K2*(CA_eq(H, params)-2*params.DIC)

function pCO₂(DIC, ALK, T, S, params)
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    ALK *= 1.e-3/params.ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC *= 1.e-3/params.ρₒ # micromoles to moles    

    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity

    K0 = exp(-60.2409 + 9345.17/T + 23.3585*log(T/100) + S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/T - 1.468591*log(T) + (-0.138681 - 9.33291/T)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/T -1.437139*log(T) + (-0.128417 - 24.41239/T)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/T + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T) + 0.053105*sqrt(S)*T)

    H = 10^(-params.pH) # initial guess from arg list

    p=(DIC=DIC, ALK=ALK, K0=K0, K1=K1, K2=K2, KB=KB, Boron=Boron)

    H = find_zero(H_eq, H, atol=1e-100, p=p)
    CA = CA_eq(H, p)

    CO2aq = CA/(K1/H + 2*K1*K2/H^2)*1e6 # Eq 11  μmol/kg
    return CO2aq/K0 #ppm
end

#####
##### Gas exchange velocity
#####

k(gas::Symbol, T, params)=0.39*(0.01/3600)*params.uₐᵥ^2*(Sc(T, getproperty(params.Sc_params, gas))/660)^(-0.5)#m/s, may want to add variable wind speed instead of average wind here at some point
Sc(T, params) = params.A-params.B*T+params.C*T^2-params.D*T^3

α(gas::Symbol, T, S, params)=β(T+273.15, S, getproperty(params.β_params, gas))*(T+273.15)*0.00367#/(T+273.15) - disagree with origional paper but this matches dimensionless Henry coefficiet of 3.2x10⁻² at 298.15K, S=0. See https://www.wikiwand.com/en/Henry%27s_law 
β(T, S, params) = exp(params.A₁+params.A₂*(100/T)+params.A₃*log(T/100)+S*(params.B₁+params.B₂*(T/100)+params.B₃*(T/100)^2))

#now fairly sure that this is giving the correct result for CO₂  as βρ is the henrys coefficient which sould be ∼34mol/m³ atm
K(gas::Symbol, T, S, params)=k(gas, T, params)*β(T+273.15, S, getproperty(params.β_params, gas))*params.ρₒ #L=ρ\_wK₀ ->  https://reader.elsevier.com/reader/sd/pii/0304420374900152 and here K₀=β

#####
##### Air sea gas flux boundary functions
#####

function airseaflux(x, y, t, T::AbstractFloat, S::AbstractFloat, conc::AbstractFloat, params)
    #https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/92JC00188 could add chemical enhancement factor
    if params.gas in (:O₂, )#doing like this for flexability in the future
        return k(params.gas, T, params)*(conc-α(params.gas, T, S, params)*getproperty(params.conc_air, params.gas))
    elseif params.gas in (:CO₂, )
        return K(params.gas, T, S, params)*(conc-getproperty(params.conc_air, params.gas))/(1000)#mol*ppmv/m^2*atm*s to mmolC/m²s not sure this is correct
    else
        throw(ArgumentError("Invalid gas choice for airseaflux"))
    end
end

function airseaflux(x, y, t, DIC::AbstractFloat, ALK::AbstractFloat, T::AbstractFloat, S::AbstractFloat, params)
    if !(params.gas in (:CO₂, ))
        throw(ArgumentError("Too many arguments for gas type $(params.gas)"))
    end
    conc = pCO₂(DIC, ALK, T+273.15, S, params)
    return airseaflux(x, y, t, T, S, conc, params)
end

airseaflux(x, y, t, conc, params) = airseaflux(x, y, t, params.T(x, y, 0, t), params.S(x, y, 0, t), conc, params)
airseaflux(x, y, t, DIC, ALK, params) = airseaflux(x, y, t, DIC, ALK, params.T(x, y, 0, t), params.S(x, y, 0, t), params)

#####
##### Boundary condition setup
#####
"""
    airseasetup(gas::Symbol; forcing=(T=nothing, S=nothing), parameters=defaults.airseaflux)

Returns an Oceananigans `FluxBoundaryCondition` for the tracer relivant to the specified gas (i.e. DIC for CO₂ and oxygen for O₂).
Arguments
=========
* `gas`: Symbol specifying the gas to exchange into the water, current choices are `:CO₂` and `:O₂`
Keyword arguments
=================
* `forcing`: NamedTuple of functions in the form `func(x, y, z, t)` for temperature, `T`, and salinity `S`
    Optional as the fallback is to look for tracer fields by these names (the prefered method where the physics are being resolved).
* `parameters`: NamedTuple of parameters to replace the defaults
"""

function airseasetup(gas::Symbol; forcings=(T=nothing, S=nothing), parameters=defaults.airseaflux)
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