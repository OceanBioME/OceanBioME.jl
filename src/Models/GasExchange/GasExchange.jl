"""
`GasExchangeModel` to solve chemical equilibrium parameterisations
"""
module GasExchangeModel

export GasExchange, 
       CarbonDioxideGasExchangeBoundaryCondition, 
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

# fallbacks
field_dependencies(name) = tuple()
optional_fields(name) = tuple()

# just get the value from the model tracer
struct TracerConcentration{T}
    tracer :: T
end

summary(tc::TracerConcentration) = "Model tracer concentration of $(tc.tracer)"
show(io::IO, tc::TracerConcentration) = println(io, summary(tc))

@inline surface_value(::TracerConcentration, i, j, grid, clock, args...) = @inbounds args[end][i, j, grid.Nz]

field_dependencies(mc::TracerConcentration) = (:T, mc.tracer)

# include all the bits
include("surface_values.jl")
include("generic_parameterisations.jl")
include("gas_exchange.jl")
include("carbon_dioxide_concentration.jl")
include("schmidt_number.jl")
include("gas_transfer_velocity.jl")

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

`water_concentration` should usually be a `TracerConcentration` (which specifies to read
the concentration directly from the model), or a `CarbonDioxideConcentration` which diagnoses
the partial pressure of CO₂ in the water.

`transfer_velocity` should be a function of the form `k(u₁₀, T)`.
"""
function GasExchangeBoundaryCondition(; water_concentration,
                                        air_concentration,
                                        transfer_velocity,
                                        wind_speed,
                                        discrete_form = false)

    wind_speed = normalise_surface_function(wind_speed; discrete_form)
    air_concentration = normalise_surface_function(air_concentration; discrete_form)

    exchange_function = GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration)

    return FluxBoundaryCondition(exchange_function; discrete_form = true)
end

"""
    CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry = CarbonChemistry(),
                                                transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = CarbonDioxidePolynomialSchmidtNumber()),
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
function CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry = CarbonChemistry(),
                                                     transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = CarbonDioxidePolynomialSchmidtNumber()),
                                                     air_concentration = 413, # ppmv
                                                     wind_speed = 2,
                                                     water_concentration = nothing,
                                                     silicate_and_phosphate_names = nothing,
                                                     kwargs...)

    if isnothing(water_concentration)
        water_concentration = CarbonDioxideConcentration(; carbon_chemistry, silicate_and_phosphate_names)
    elseif !isnothing(carbon_chemistry)
        @warn "Make sure that the `carbon_chemistry` $(carbon_chemistry) is the same as that in `water_concentration` $(water_concentration) (or set it to `nothing`)"
    end

    return GasExchangeBoundaryCondition(; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)
end

"""
    OxygenGasExchangeBoundaryCondition(; transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = OxygenPolynomialSchmidtNumber()),
                                         water_concentration = TracerConcentration(:O₂),
                                         air_concentration = 9352.7, # mmolO₂/m³
                                         wind_speed = 2,
                                         kwagrs...)

Returns a `FluxBoundaryCondition` for the gas exchange between oxygen dissolved in the water
specified by the the `TracerConcentration` in the base model, and `air_concentration` with `transfer_velocity`
(see `GasExchangeBoundaryCondition` for details).

`kwargs` are passed on to `GasExchangeBoundaryCondition`.
"""
OxygenGasExchangeBoundaryCondition(; transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = OxygenPolynomialSchmidtNumber()),
                                     water_concentration = TracerConcentration(:O₂),
                                     air_concentration = 9352.7, # mmolO₂/m³
                                     wind_speed = 2,
                                     kwargs...) = 
    GasExchangeBoundaryCondition(; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)

end # module