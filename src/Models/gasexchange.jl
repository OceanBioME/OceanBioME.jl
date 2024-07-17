#####
##### Gas exchange model of [Wanninkhof1992](@citet)
#####
# TODO: Implement Ho et al. 2006 wind speed dependence

k(T, uₐᵥ, Sc_params) = 0.39 * (0.01 / 3600) * uₐᵥ ^ 2 * (Sc(T, Sc_params) / 660) ^ (-0.5)# m/s, may want to add variable wind speed instead of average wind here at some point
Sc(T, params) = params.A - params.B * T + params.C * T ^ 2 - params.D * T ^ 3

α(T, S, β_params) = β(T + 273.15, S, β_params) * (T + 273.15) * 0.00367#/(T+273.15) - disagree with origional paper but this matches dimensionless Henry coefficiet of 3.2x10⁻² at 298.15K, S=0. See https://www.wikiwand.com/en/Henry%27s_law 
β(T, S, params) = exp(params.A₁ + params.A₂ * (100 / T) + params.A₃ * log(T / 100) + S * (params.B₁ + params.B₂ * (T / 100) + params.B₃ * (T / 100) ^ 2))

#now fairly sure that this is giving the correct result for CO₂  as βρ is the henrys coefficient which sould be ∼34mol/m³ atm
K(T, S, uₐᵥ, Sc_params, β_params, ρₒ) = k(T, uₐᵥ, Sc_params) * β(T + 273.15, S, β_params) * ρₒ #L=ρ\_wK₀ ->  https://reader.elsevier.com/reader/sd/pii/0304420374900152 and here K₀=β

#####
##### Boundary condition setup
#####

struct GasExchange{G, ScP, βP, FT, AC, AP, PCO} <: Function
    gas :: G

    schmidt_params :: ScP
    solubility_params :: βP
    ocean_density :: FT
    air_concentration :: AC
    air_pressure :: AP
    average_wind_speed :: FT

    pCO₂ :: PCO
end

adapt_structure(to, gasexchange::GasExchange) = GasExchange(adapt(to, gasexchange.gas),
                                                            adapt(to, gasexchange.schmidt_params),
                                                            adapt(to, gasexchange.solubility_params),
                                                            gasexchange.ocean_density,
                                                            adapt(to, gasexchange.air_concentration),
                                                            adapt(to, gasexchange.air_pressure),
                                                            gasexchange.average_wind_speed,
                                                            adapt(to, gasexchange.pCO₂))

"""
    GasExchange(; gas,
                  schmidt_params::ScP = (CO₂ = (A=2073.1, B=125.62, C=3.6276, D=0.043219),
                                   O₂ = (A=1953.4, B=128.0, C=3.9918, D=0.050091))[gas],
                  solubility_params::βP = (CO₂ = (A₁=-60.2409, A₂=93.4517, A₃=23.3585, B₁=0.023517, B₂=-0.023656, B₃=0.0047036),
                                     O₂ = (A₁=-58.3877, A₂=85.8079, A₃=23.8439, B₁=-0.034892, B₂=0.015568, B₃=-0.0019387))[gas],
                  ocean_density::FT = 1026, # kg/m³
                  air_concentration::AC = (CO₂ = 413.4, O₂ = 9352.7)[gas], # ppmv, mmolO₂/m³ (20.95 mol O₂/mol air, 0.0224m^3/mol air)
                  air_pressure::FT = 1.0, # atm
                  average_wind_speed::FT = 10, # m/s
                  field_dependencies = (CO₂ = (:DIC, :ALK), O₂ = (:OXY, ))[gas])

Construct an Oceananigans `FluxBoundaryCondition` for the exchange of `gas` with the relevant tracer (i.e., DIC for CO₂ and oxygen for O₂).
Please see note for other gases.

Keyword arguments
=================

- `gas`: (required) the gas to be exchanged, if `:CO₂` or `:O₂` are specified then all other settings may be infered
- `schmidt_params` : named tuple of parameters for calculating the Schmidt number using the parameterisation of [Wanninkhof1992](@citet)
- `solubility_params` : named tuple of parameters for calculating the solubility (for O₂ the Bunsen solubility and CO₂ K₀, see note)
- `ocean_density` : density of the ocean in kg/m³
- `air_concentratio` : concentration of the gas in air in relivant units (i.e. ppmv for CO₂ and mmol O₂/m³ for O₂), can also be a function of x, y, t, or a field
- `air_pressure` : air pressure in atm (only used for CO₂), can also be a function of x, y, t, or a field
- `average_wind_speed` : average wind speed at 10m used to calculate the gas transfer velocity by the [Wanninkhof1992](@citet) parameterisation
- `field_dependencies` : tracer fields that gas exchange depends on, if the defaults have different names in your model you can specify as long as they are in the same order
- `pCO₂` : pCO₂ calculator

!!! note "Gases other than CO₂ and O₂"
    This model is fully capable of exchanging any gas but the parameters have only been configured for CO₂ and O₂, and the specific formulation
    is only ensured for these gasses. For any gas where the [Wanninkhof1992](@citet) parameterisation returns the Bunsen Solubility Coefficient
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
                       pCO₂::PCO = gas == :CO₂ ? CarbonChemistry() : nothing) where {ScP, βP, FT, AC, AP, PCO}

    gas = Val(gas)
    G = typeof(gas)

    gasexchange =  GasExchange{G, ScP, βP, FT, AC, AP, PCO}(gas, 
                                                            schmidt_params, 
                                                            solubility_params, 
                                                            ocean_density, 
                                                            air_concentration, 
                                                            air_pressure, 
                                                            average_wind_speed, 
                                                            pCO₂)

    return FluxBoundaryCondition(gasexchange, field_dependencies = (field_dependencies..., :T, :S))
end

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
#=It is not possible for this to work on GPU since we need the grid and locs before (since fields are reduced to vectors)
  We probably need a more generic and better way todo this but here is not the place.
  For now if you wanted to use data you could do similar to https://github.com/OceanBioME/GlobalOceanBioME.jl/blob/main/src/one_degree_surface_par.jl
# interpolate doesn't really work on 2D fields
@inline function get_value(x, y, t, conc::Field{LX, LY, LZ}) where {LX, LY, LZ}
    grid = conc.grid

    i, j, _ = fractional_indices((x, y, 0.0), grid, LX(), LY(), Center())

    ξ, i = mod(i, 1), Base.unsafe_trunc(Int, i)
    η, j = mod(j, 1), Base.unsafe_trunc(Int, j)

    _, _, ks = indices(conc)

    return @inbounds  ((1 - ξ) * (1 - η) * conc[i,   j,   ks[1]]
                     + (1 - ξ) *      η  * conc[i,   j+1, ks[1]]
                     +      ξ  * (1 - η) * conc[i+1, j,   ks[1]]
                     +      ξ  *      η  * conc[i+1, j+1, ks[1]])
end
=#