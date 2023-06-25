##### 
##### Carbonate chemistry to determine pCO₂
##### As per OCMIP Phase 2 http://ocmip5.ipsl.jussieu.fr/OCMIP/phase2/simulations/Abiotic/Cchem/co2calc.f simplified as in PISCESv2
#####

struct pCO₂{P0, P1, P2, PB, PW, FT}
                  solubility :: P0
    bicarbonate_dissociation :: P1
      carbonate_dissociation :: P2
     boric_acid_dissociation :: PB 
          water_dissociaiton :: PW # don't think this is the right name for it
              lower_pH_bound :: FT
              upper_pH_bound :: FT
                 boron_ratio :: FT
           thermal_expansion :: FT
          haline_contraction :: FT


    pCO₂(solubility::P0,
         bicarbonate_dissociation::P1,
         carbonate_dissociation::P2,
         boric_acid_dissociation::PB,
         water_dissociaiton::PW,
         lower_pH_bound::FT,
         upper_pH_bound::FT,
         boron_ratio::FT,
         thermal_expansion::FT,
         haline_contraction::FT) where {P0, P1, P2, PB, PW, FT} =
         new{P0, P1, P2, PB, PW, FT}(solubility, 
                                     bicarbonate_dissociation, 
                                     carbonate_dissociation, 
                                     boric_acid_dissociation, 
                                     water_dissociaiton, 
                                     lower_pH_bound, upper_pH_bound, 
                                     boron_ratio, 
                                     thermal_expansion, haline_contraction)
end

adapt_structure(to, pCO₂_model::pCO₂) = pCO₂(adapt(to, pCO₂_model.solubility),
                                             adapt(to, pCO₂_model.bicarbonate_dissociation),
                                             adapt(to, pCO₂_model.carbonate_dissociation),
                                             adapt(to, pCO₂_model.boric_acid_dissociation),
                                             adapt(to, pCO₂_model.water_dissociaiton),
                                             pCO₂_model.lower_pH_bound, pCO₂_model.upper_pH_bound,
                                             pCO₂_model.boron_ratio, pCO₂_model.thermal_expansion, pCO₂_model.haline_contraction)

@inline function titrate_alkalinity(H, p)
    return p.DIC * (p.k¹ * H + 2 * p.k¹ * p.k²) / (H ^ 2 + p.k¹ * H + p.k¹ * p.k²) -
           (p.Alk - p.kʷ / H - p.boron / (1 + H / p.kᵇ))
end

@inline function (p::pCO₂)(DIC, Alk, T, S)
    ρₒ = 1027 * (1 - p.thermal_expansion * T + p.haline_contraction * S)

    T += 273.15

    Alk *= 1e-3 / ρₒ
    DIC *= 1e-3 / ρₒ

    pk⁰ = p.solubility
    pk¹ = p.bicarbonate_dissociation
    pk² = p.carbonate_dissociation
    pkᵇ = p.boric_acid_dissociation
    pkʷ = p.water_dissociaiton

    k¹ = 10 ^ (pk¹.C + pk¹.invT / T + pk¹.logT * log(T) + pk¹.S * S + pk¹.S² * S ^ 2)

    k² = 10 ^ (pk².C + pk².invT / T + pk².S * S + pk².S² * S ^ 2)

    kᵇ = exp(pkᵇ.C + (pkᵇ.invT + pkᵇ.invTsqrtS * sqrt(S) + pkᵇ.invTS * S + pkᵇ.invTS¹⁵ * S ^ 1.5 + pkᵇ.invTS² * S ^ 2) / T
    		+ pkᵇ.sqrtS * sqrt(S) + pkᵇ.S * S + (pkᵇ.logT + pkᵇ.logTsqrtS * sqrt(S) + pkᵇ.logTS * S) * log(T) + pkᵇ.TsqrtS * sqrt(S) * T)

    kʷ = exp(pkʷ.C + pkʷ.invT / T + pkʷ.logT * log(T) +
    		 (pkʷ.sqrtS + pkʷ.sqrtSinvT / T + pkʷ.sqrtSlogT * log(T)) * sqrt(S)
             + pkʷ.S * S)

    boron = p.boron_ratio * S
    
    params = (; k¹, k², kᵇ, kʷ, DIC, Alk, boron)

    if titrate_alkalinity(10 ^ - p.upper_pH_bound, params) * titrate_alkalinity(10 ^ - p.lower_pH_bound, params) < 0
        H = find_zero(titrate_alkalinity, (10 ^ - p.upper_pH_bound, 10 ^ - p.lower_pH_bound), Bisection(); atol = 1e-10, p = params)
    else
        H = NaN
    end

    ff = exp(pk⁰.C + pk⁰.invT / T  +
             pk⁰.ClogT * log(T * pk⁰.logCT) + pk⁰.T² * T ^ 2 +
    		 S * (pk⁰.S + pk⁰.ST * T + pk⁰.ST² * T ^ 2))

    CO₂ = DIC * H ^ 2/ (H ^ 2 + k¹ * H + k¹ * k²)
    pCO₂ = (CO₂ / ff) * 10 ^ 6

    return pCO₂ # μatm
end

OCMIP_solubility = (C = -162.8301, invT = 218.2968 * 100, logCT = 1 / 100, ClogT = 90.9241, T² = - 1.47696 / (100 ^ 2), ST² = 0.0049867 / (100 ^ 2), ST = -0.025225 / 100, S = .025695)#

OCMIP_bicarbonate_dissociation = (C = 62.008, S = 0.0118, S² = -0.000116, invT = -3670.7, logT = -9.7944)

OCMIP_carbonate_dissociation = (C = -4.777, S = 0.0184, S² = -0.000118, invT = -1394.7, logT = 0.0)

OCMIP_boric_acid_dissociation = (C = 148.0248, invT = -8966.9, invTsqrtS = -2890.53, invTS = -77.942, invTS¹⁵ = 1.728, invTS² = - 0.0996,
                                 sqrtS = 137.1942, S = 1.62142, logT = - 24.4344, logTsqrtS = - 25.085, logTS = - 0.2474, TsqrtS = 0.053105)

OCMIP_water_dissociaiton = (C = 148.9652, invT = - 13847.26, logT = - 23.6521, sqrtSinvT = 118.67, sqrtS = -5.977, sqrtSlogT = 1.0495, S = -0.01615)

OCMIP_default = pCO₂(OCMIP_solubility, 
                     OCMIP_bicarbonate_dissociation, 
                     OCMIP_carbonate_dissociation, 
                     OCMIP_boric_acid_dissociation, 
                     OCMIP_water_dissociaiton, 
                     0.0, 14.0, 
                     0.000232 / 1.80655 / 10.811,
                     1.67e-4, 7.80e-4)

#####
##### Gas exchange model of [Wanninkhof1992](@cite)
#####
# TODO: Impliment Ho et al. 2006 wind speed dependence

k(T, uₐᵥ, Sc_params) = 0.39 * (0.01 / 3600) * uₐᵥ ^ 2 * (Sc(T, Sc_params) / 660) ^ (-0.5)# m/s, may want to add variable wind speed instead of average wind here at some point
Sc(T, params) = params.A - params.B * T + params.C * T ^ 2 - params.D * T ^ 3

α(T, S, β_params) = β(T + 273.15, S, β_params) * (T + 273.15) * 0.00367#/(T+273.15) - disagree with origional paper but this matches dimensionless Henry coefficiet of 3.2x10⁻² at 298.15K, S=0. See https://www.wikiwand.com/en/Henry%27s_law 
β(T, S, params) = exp(params.A₁ + params.A₂ * (100 / T) + params.A₃ * log(T / 100) + S * (params.B₁ + params.B₂ * (T / 100) + params.B₃ * (T / 100) ^ 2))

#now fairly sure that this is giving the correct result for CO₂  as βρ is the henrys coefficient which sould be ∼34mol/m³ atm
K(T, S, uₐᵥ, Sc_params, β_params, ρₒ) = k(T, uₐᵥ, Sc_params) * β(T + 273.15, S, β_params) * ρₒ #L=ρ\_wK₀ ->  https://reader.elsevier.com/reader/sd/pii/0304420374900152 and here K₀=β

#####
##### Boundary condition setup
#####

struct GasExchange{G, ScP, βP, FT, AC, AP, T, S, PCO}
    gas :: G

    schmidt_params :: ScP
    solubility_params :: βP
    ocean_density :: FT
    air_concentration :: AC
    air_pressure :: AP
    average_wind_speed :: FT

    temperature :: T
    salinity :: S

    pCO₂ :: PCO
end

adapt_structure(to, gasexchange::GasExchange) = GasExchange(adapt(to, gasexchange.gas),
                                                            adapt(to, gasexchange.schmidt_params),
                                                            adapt(to, gasexchange.solubility_params),
                                                            gasexchange.ocean_density,
                                                            adapt(to, gasexchange.air_concentration),
                                                            adapt(to, gasexchange.air_pressure),
                                                            gasexchange.average_wind_speed,
                                                            adapt(to, gasexchange.temperature),
                                                            adapt(to, gasexchange.salinity),
                                                            adapt(to, gasexchange.pCO₂))

"""
    GasExchange(;gas,
                schmidt_params::ScP = (CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219),
                                O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091))[gas],
                solubility_params::βP = (CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036),
                                    O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387))[gas],
                ocean_density::FT = 1026, # kg/m³
                air_concentration::AC = (CO₂ = 413.4, O₂ = 9352.7)[gas], # ppmv, mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
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
- `ocean_density` : density of the ocean in kg/m³
- `air_concentratio` : concentration of the gas in air in relivant units (i.e. ppmv for CO₂ and mmol O₂/m³ for O₂), can also be a function of x, y, t, or a field
- `air_pressure` : air pressure in atm (only used for CO₂), can also be a function of x, y, t, or a field
- `average_wind_speed` : average wind speed at 10m used to calculate the gas transfer velocity by the [Wanninkhof1992](@cite) parameterisation
- `field_dependencies` : tracer fields that gas exchange depends on, if the defaults have different names in your model you can specify as long as they are in the same order
- `temperature` : either `nothing` to track a temperature tracer field, or a function or shape `f(x, y, z, t)` for the temperature in °C
- `salinity` : either `nothing` to track a salinity tracer field, or a function or shape `f(x, y, z, t)` for the salinity in ‰
- `pCO₂` : pCO₂ calculator

Note
====

This model is fully capable of exchanging any gas but the parameters have only been configured for CO₂ and O₂, and the specific formulaiton
is only ensured for these gasses. For any gas where the [Wanninkhof1992](@cite) parameterisation returns the Bunsen Solubility Coefficient
this model will work out of the box and can just be passed new parameters. For the other solubility types (i.e. K₀, K' and f) you will need
to overload the `(gasexchange::GasExchange)` function to ensure the correct formulaiton.
"""

function GasExchange(; gas,
                       schmidt_params::ScP = (CO₂ = (A = 2073.1, B = 125.62, C = 3.6276, D = 0.043219),
                                              O₂ = (A = 1953.4, B = 128.0, C = 3.9918, D = 0.050091))[gas],
                       solubility_params::βP = (CO₂ = (A₁ = -60.2409, A₂ = 93.4517, A₃ = 23.3585, B₁ = 0.023517, B₂ = -0.023656, B₃ = 0.0047036),
                                                O₂ = (A₁ = -58.3877, A₂ = 85.8079, A₃ = 23.8439, B₁ = -0.034892, B₂ = 0.015568, B₃ = -0.0019387))[gas],
                       ocean_density::FT = 1024.5, # kg/m³
                       air_concentration::AC = (CO₂ = 413.4, O₂ = 9352.7)[gas], # ppmv, mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
                       air_pressure::AP = 1.0, # atm
                       average_wind_speed::FT = 10.0, # m/s
                       field_dependencies = (CO₂ = (:DIC, :Alk), O₂ = (:O₂, ))[gas],
                       temperature::T = nothing,
                       salinity::S = nothing,
                       pCO₂::PCO = gas == :CO₂ ? OCMIP_default : nothing) where {ScP, βP, FT, AC, AP, T, S, PCO}

    gas = Val(gas)
    G = typeof(gas)

    gasexchange =  GasExchange{G, ScP, βP, FT, AC, AP, T, S, PCO}(gas, 
                                                                  schmidt_params, 
                                                                  solubility_params, 
                                                                  ocean_density, 
                                                                  air_concentration, 
                                                                  air_pressure, 
                                                                  average_wind_speed, 
                                                                  temperature, 
                                                                  salinity,
                                                                  pCO₂)

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
    conc = gasexchange.pCO₂(DIC, ALK, T, S)
    return gasexchange(x, y, t, conc, T, S)
end

@inline function (gasexchange::GasExchange)(x, y, t, conc, T, S) 
    return k(T, gasexchange.average_wind_speed, gasexchange.schmidt_params) * (conc - α(T, S, gasexchange.solubility_params) * get_value(x, y, t, gasexchange.air_concentration))
end

@inline function (gasexchange::GasExchange{<:Val{:CO₂}, <:Any, <:Any, <:Any, <:Any, <:Any})(x, y, t, conc, T, S) 
    return K(T, S, gasexchange.average_wind_speed, gasexchange.schmidt_params, gasexchange.solubility_params, gasexchange.ocean_density) * (conc - get_value(x, y, t, gasexchange.air_concentration) * get_value(x, y, t, gasexchange.air_pressure)) / 1000#μmol/m²s to mmolC/m²s not sure this is correct
end

@inline get_value(x, y, t, air_concentration::Number) = air_concentration
@inline get_value(x, y, t, air_concentration::Function) = air_concentration(x, y, t)

# interpolation does not work on 2d grids as trilinear interpolation fails but we're only ever going to be asking for values on grid points
@inline function get_value(x, y, t, air_concentration::Field{Center, Center, Center})
    (ξ, i), (η, j), (ζ, k) = modf.(fractional_indices(x, y, znodes(air_concentration)[end], (Center(), Center(), Center()), air_concentration.grid))
    return air_concentration[floor(Int, i) + 1, floor(Int, j) + 1, floor(Int, k) + 1]
end