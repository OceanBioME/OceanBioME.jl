module AirSeaFlux
using Roots
function dic(x, y, t, DIC, ALK, T::AbstractFloat, S::AbstractFloat, params) # has to be only including x,y,t without z, because this will apply to the z direction. f(x, y, t) on z-boundaries.
    #https://clima.github.io/OceananigansDocumentation/stable/model_setup/boundary_conditions/
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    #ALK *= 1.e-6 # microequivalents to equivalents
    #DIC *= 1.e-6 # micromoles to moles
    ALK *= 1.e-3/params.ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC *= 1.e-3/params.ρₒ # micromoles to moles    

    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity
    
    K0 = exp(-60.2409 + 9345.17/T + 23.3585*log(T/100) + S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/T - 1.468591*log(T) + (-0.138681 - 9.33291/T)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/T -1.437139*log(T) + (-0.128417 - 24.41239/T)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/T + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T) + 0.053105*sqrt(S)*T)

    H = 10^(-params.pH) # initial guess from arg list
    #potential peformance botteneck
    @inline CA(H) = ALK - (KB/(KB + H))*Boron
    @inline H_eq(H) = CA(H)*H^2 + K1*(CA(H)-DIC)*H + K1*K2*(CA(H)-2*DIC)
    H = find_zero(H_eq, H, atol=1.e-15)
    
    #pH = -log10(H)
    CO2aq = CA(H)/(K1/H + 2*K1*K2/H^2)*1e6 # Eq 11  μmol/kg
    pCO2 = CO2aq/K0 # Eq 4 (converted from atm to ppmv) ppm 
    HCO3 = K0*K1*CO2aq/K0/H       #μmol/kg
    CO3 = K0*K1*K2*CO2aq/K0/H^2   #μmol/kg
    DIC = CO2aq + HCO3 + CO3
    R = DIC/CO3
    #[flux, pCO2, CO2aq, HCO3, CO3, DIC, R]
    flux = 7.7e-4*params.U_10^2*(pCO2-params.pCO2_air)/(365*24*3600/1000) # mmol/m^2/S

    return flux      #Wanninkhof 2014 equ.6 positive value means upward flux meaning losing Carbon 
end

dic(x, y, t, DIC, ALK, params) = dic(x, y, t, DIC, ALK, params.T(x, y, 0, t), params.S(x, y, 0, t), params)

#not the nicest or most julian way todo this but don't know how else to destinguish them
function dic(x, y, t, DIC, ALK, TorS, params)
    if :T in keys(params)
        return dic(x, y, t, DIC, ALK, params.T(x, y, 0, t), TorS, params)
    else 
        return dic(x, y, t, DIC, ALK, TorS, params.S(x, y, 0, t), params)
    end
end
end