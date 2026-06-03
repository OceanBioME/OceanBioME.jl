"""
`GasExchangeModel` to solve chemical equilibrium parameterisations
"""
module GasExchangeModel

export GasExchange, 
       CarbonDioxideGasExchangeBoundaryCondition,
       CarbonDioxideGasExchangeBoundaryConditions, 
       OxygenGasExchangeBoundaryCondition, 
       GasExchangeBoundaryCondition,
       ScaledGasTransferVelocity,
       SchmidtScaledTransferVelocity,
       CarbonDioxidePolynomialSchmidtNumber,
       OxygenPolynomialSchmidtNumber

using Adapt
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Center
using Oceananigans.Grids: xnode, ynode

using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

import Base: show, summary
import Adapt: adapt_structure

const ATM = 101325 # Pa
const GAS_CONSTANT = 8.31446261815324 # J / kg / mol

include("surface_values.jl")

# just get the value from the model tracer
struct OxygenConcentration end

summary(::OxygenConcentration) = "Model tracer `OxygenConcentration`"

@inline surface_value(::OxygenConcentration, i, j, grid, clock, model_fields) = @inbounds model_fields.O₂[i, j, grid.Nz]

# include all the bits
include("generic_parameterisations.jl")
include("gas_exchange.jl")
include("carbon_dioxide_concentration.jl")
include("schmidt_number.jl")
include("gas_transfer_velocity.jl")
include("gas_solubility.jl")

using .ScaledGasTransferVelocity

# wrappers to produce boundary conditions

"""
    GasExchangeBoundaryCondition(; water_concentration,
                                   air_concentration,
                                   transfer_velocity,
                                   wind_speed)

Returns a `FluxBoundaryCondition` for the gas exchange between `water_concentration` and `air_concentration`
with `transfer_velocity`.

`water_concentration`, `air_concentration` and `wind_speed` can either be numbers, 
functions of the form `(x, y, t)`, functions of the form `(i, j, grid, clock, model_fields)` 
if `discrete_form` is set to true, or any kind of `Field`.

`water_concentration` should usually be a `[Tracer]Concentration` where is the name of the
tracer (you will have to build your own if this is not `OxygenConcentration`), 
or a `CarbonDioxideConcentration` which diagnoses the partial pressure of CO₂ in the water.

`transfer_velocity` should be a function of the form `k(u₁₀, T)`.
"""
function GasExchangeBoundaryCondition(FT = Float64; 
                                      water_concentration,
                                      air_concentration,
                                      transfer_velocity,
                                      wind_speed,
                                      discrete_form = false)

    wind_speed = normalise_surface_function(wind_speed; discrete_form, FT)
    air_concentration = normalise_surface_function(air_concentration; discrete_form, FT)

    exchange_function = GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration)

    return FluxBoundaryCondition(exchange_function; discrete_form = true)
end

"""
    CarbonDioxideGasExchangeBoundaryCondition(FT = Float64; 
                                              carbon_chemistry = CarbonChemistry(FT),
                                              transfer_velocity = SchmidtScaledTransferVelocity(schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT)),
                                              air_concentration = 413, # ppmv
                                              wind_speed = 2,
                                              water_concentration = nothing,
                                              silicate_and_phosphate_names = nothing,
                                              kwargs...)

Returns a `FluxBoundaryCondition` for the gas exchange between carbon dioxide dissolved in the water
specified by the `carbon_chemisty` model, and `air_concentration` with `transfer_velocity` (see 
`GasExchangeBoundaryCondition` for details).

`silicate_and_phosphate_names` should either be `nothing`, a `Tuple`` of symbols specifying the name of the silicate
and phosphate tracers, or a `NamedTuple`  of values for the `carbon_chemistry` model.

`kwargs` are passed on to `GasExchangeBoundaryCondition`.

Note: The model always requires `T`, `S`, `DIC`, and `Alk` to be present in the model.
"""
CarbonDioxideGasExchangeBoundaryCondition(FT = Float64; 
                                          carbon_chemistry = CarbonChemistry(FT),
                                          transfer_velocity = 
                                               SchmidtScaledTransferVelocity(FT; 
                                                  schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT),
                                                  solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(carbon_chemistry.solubility,
                                                                                                        carbon_chemistry.density_function)),
                                          air_concentration = 413, # ppmv
                                          wind_speed = 2,
                                          water_concentration = nothing,
                                          silicate_and_phosphate_names = nothing,
                                          kwargs...) =
    CarbonDioxideGasExchangeBoundaryConditions(1, FT;
                                               carbon_chemistry,
                                               transfer_velocity,
                                               air_concentration,
                                               wind_speed,
                                               water_concentration,
                                               silicate_and_phosphate_names,
                                               kwargs...)


"""
    CarbonDioxideGasExchangeBoundaryConditions(N = 1, FT = Float64; 
                                               carbon_chemistry = CarbonChemistry(FT),
                                               transfer_velocity = SchmidtScaledTransferVelocity(schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT)),
                                               air_concentration = 413, # ppmv
                                               wind_speed = 2,
                                               water_concentration = nothing,
                                               silicate_and_phosphate_names = nothing,
                                               kwargs...)

Returns a `FluxBoundaryCondition` for the gas exchange between carbon dioxide dissolved in the water
specified by the `carbon_chemisty` model, and `air_concentration` with `transfer_velocity` (see 
`GasExchangeBoundaryCondition` for details).

`silicate_and_phosphate_names` should either be `nothing`, a `Tuple`` of symbols specifying the name of the silicate
and phosphate tracers, or a `NamedTuple`  of values for the `carbon_chemistry` model.

`kwargs` are passed on to `GasExchangeBoundaryCondition`.

Note: The model always requires `T`, `S`, `DIC`, and `Alk` to be present in the model.
"""
function CarbonDioxideGasExchangeBoundaryConditions(N = 1, FT = Float64; 
                                                    carbon_chemistry = CarbonChemistry(FT),
                                                    transfer_velocity = 
                                                         SchmidtScaledTransferVelocity(FT; 
                                                            schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT),
                                                            solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(carbon_chemistry.solubility,
                                                                                                                  carbon_chemistry.density_function)),
                                                    air_concentration = 413, # ppmv
                                                    wind_speed = 2,
                                                    water_concentration = nothing,
                                                    silicate_and_phosphate_names = nothing,
                                                    kwargs...)

    if isnothing(water_concentration)
        if N == 1
            water_concentration = CarbonDioxideConcentration(FT; carbon_chemistry, silicate_and_phosphate_names)
        else
            water_concentration = map(
                n->CarbonDioxideConcentration(FT; carbon_chemistry, 
                                                  silicate_and_phosphate_names,
                                                  DIC = Symbol(:DIC, n),
                                                  Alk = Symbol(:Alk, n)),
                1:N)
        end
    end

    if N == 1
        return GasExchangeBoundaryCondition(FT; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)
    else
        dic_names = tuple(map(n->Symbol(:DIC, n), 1:N)...)
        dic_boundary_kwargs = (; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)
        return NamedTuple{dic_names}(map(n->materialise_multiple_dic_boundaries(n, FT; dic_boundary_kwargs...), 1:N))
    end
end

materialise_multiple_dic_boundaries(n, FT; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...) =
    GasExchangeBoundaryCondition(FT; water_concentration = one_or_nth(water_concentration, n),
                                     air_concentration = one_or_nth(air_concentration, n),
                                     transfer_velocity = one_or_nth(transfer_velocity, n),
                                     wind_speed = one_or_nth(wind_speed, n),
                                     kwargs...)

one_or_nth(only_one, n) = only_one
one_or_nth(n_options::Array, n) = n_options[n]
one_or_nth(n_options::Tuple, n) = n_options[n]

"""
    OxygenGasExchangeBoundaryCondition(FT = Float64; 
                                       transfer_velocity = SchmidtScaledTransferVelocity(schmidt_number = OxygenPolynomialSchmidtNumber(FT)),
                                       water_concentration = OxygenConcentration(),
                                       air_concentration = 9352.7, # mmolO₂/m³
                                       wind_speed = 2,
                                       kwagrs...)

Returns a `FluxBoundaryCondition` for the gas exchange between oxygen dissolved in the water
specified by the the `OxygenConcentration` in the base model, and `air_concentration` with `transfer_velocity`
(see `GasExchangeBoundaryCondition` for details).

`kwargs` are passed on to `GasExchangeBoundaryCondition`.
"""
OxygenGasExchangeBoundaryCondition(FT = Float64;
                                   transfer_velocity = SchmidtScaledTransferVelocity(FT; schmidt_number = OxygenPolynomialSchmidtNumber(FT)),
                                   water_concentration = OxygenConcentration(),
                                   air_concentration = PartiallySolubleGas(FT; air_concentration = 9352.7, solubility = OxygenSolubility(FT)), # mmolO₂/m³
                                   wind_speed = 2,
                                   kwargs...) = 
    GasExchangeBoundaryCondition(FT; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)

end # module
