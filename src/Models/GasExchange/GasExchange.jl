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
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

import Base: show, summary
import Adapt: adapt_structure

const ATM = 101325 # Pa
const GAS_CONSTANT = 8.31446261815324 # J / kg / mol

# fallbacks
field_dependencies(name) = tuple()
optional_fields(name) = tuple()

# extract the surface values
@inline surface_value(f, x, y, t, args...) = f # fallback
@inline surface_value(f::Function, x, y, t, args...) = f(x, y, t, args...)

struct TracerConcentration{T} <: Function
    tracer :: T
end

summary(tc::TracerConcentration) = "Model tracer concentration of $(tc.tracer)"
show(io::IO, tc::TracerConcentration) = println(io, summar(tc))

@inline (mc::TracerConcentration)(x, y, t, args...) = args[end]

field_dependencies(mc::TracerConcentration) = (:T, mc.tracer)

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

`water_concentrations` should be a type with method `water_concentrations(x, y, t, T, field_dependencies(water_concentration)...)`
(i.e. a method has to be created for `field_dependencies`).

`air_concentration` and `wind_speed` should either be numbers or functions of the form `f(x, y, t)`.

`transfer_velocity` should be a function of the form `k(u₁₀, T)`.
"""
function GasExchangeBoundaryCondition(; water_concentration,
                                        air_concentration,
                                        transfer_velocity,
                                        wind_speed)

    fd = field_dependencies(water_concentration)

    exchange_function = GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration, fd)

    return FluxBoundaryCondition(exchange_function; field_dependencies = fd)
end

"""
    CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry = CarbonChemistry(),
                                                transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = CarbonDioxidePolynomialSchmidtNumber()),
                                                air_concentration = 413, # ppmv
                                                wind_speed = 2,
                                                water_concentration = nothing,
                                                silicate_and_phosphate_names = nothing)

Returns a `FluxBoundaryCondition` for the gas exchange between carbon dioxide dissolved in the water
specified by the `carbon_chemisty` model, and `air_concentration` with `transfer_velocity`.

`silicate_and_phosphate_names` should either be `nothing` or a tuple of symbols specifying the name of the silicate
and phosphate tracers for the `carbon_chemistry` model.

Note: The model always requires `T`, `S`, `DIC`, and `Alk` to be present in the model.
"""
function CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry = CarbonChemistry(),
                                                     transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = CarbonDioxidePolynomialSchmidtNumber()),
                                                     air_concentration = 413, # ppmv
                                                     wind_speed = 2,
                                                     water_concentration = nothing,
                                                     silicate_and_phosphate_names = nothing)

    if isnothing(water_concentration)
        water_concentration = CarbonDioxideConcentration(; carbon_chemistry)
    elseif !isnothing(carbon_chemistry)
        @warn "Make sure that the `carbon_chemistry`` $(carbon_chemistry) is the same as that in `water_concentration` $(water_concentration)"
    end

    fd = (field_dependencies(water_concentration)..., ifelse(!isnothing(silicate_and_phosphate_names), silicate_and_phosphate_names, tuple())...)

    exchange_function = GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration, fd)

    return FluxBoundaryCondition(exchange_function; field_dependencies = fd)
end

"""
    OxygenGasExchangeBoundaryCondition(; transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = OxygenPolynomialSchmidtNumber()),
                                         water_concentration = TracerConcentration(:O₂),
                                         air_concentration = 9352.7, # mmolO₂/m³
                                         wind_speed = 2)

Returns a `FluxBoundaryCondition` for the gas exchange between oxygen dissolved in the water
specified by the the `TracerConcentration` in the base model, and `air_concentration` with `transfer_velocity`.
"""
OxygenGasExchangeBoundaryCondition(; transfer_velocity = SchmidtScaledTransferVelocity(; schmidt_number = OxygenPolynomialSchmidtNumber()),
                                     water_concentration = TracerConcentration(:O₂),
                                     air_concentration = 9352.7, # mmolO₂/m³
                                     wind_speed = 2) = 
    FluxBoundaryCondition(GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration, field_dependencies(water_concentration));
                          field_dependencies = (:O₂, :T))

end # module