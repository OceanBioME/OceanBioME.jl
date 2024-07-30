"""
`GasExchangeModel` to solve chemical equilibrium parameterisations
"""
module GasExchangeModel

export GasExchange, CarbonDioxideGasExchangeBoundaryCondition, OxygenGasExchangeBoundaryCondition

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

struct ModelConcentration{T} <: Function
    tracer :: T
end

@inline (mc::ModelConcentration)(x, y, t, args...) = args[end]

field_dependencies(mc::ModelConcentration) = (:T, mc.tracer)

include("gas_exchange.jl")
include("carbon_dioxide_concentration.jl")
include("schmidt_number.jl")
include("gas_transfer_velocity.jl")

# convenience constructors
function CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry = CarbonChemistry(),
                                                     transfer_velocity = ScaledTransferVelocity(; schmidt_number = CarbonDioxidePolynomialSchmidtNumber()),
                                                     air_concentration = 413, # ppmv
                                                     wind_speed = 2,
                                                     water_concentration = nothing,
                                                     include_silicate_and_carbonate = false)

    if isnothing(water_concentration)
        water_concentration = CarbonDioxideConcentration(; carbon_chemistry)
    elseif !isnothing(carbon_chemistry)
        @warn "Make sure that the `carbon_chemistry`` $(carbon_chemistry) is the same as that in `water_concentration` $(water_concentration)"
    end

    fd = (field_dependencies(water_concentration)..., ifelse(include_silicate_and_carbonate, (:silicate, :phosphate), tuple())...)
    fd_map = NamedTuple{fd}(Int[n for n in 1:length(fd)])

    exchange_function = GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration, fd, fd_map)

    return FluxBoundaryCondition(exchange_function; field_dependencies = fd)
end

OxygenGasExchangeBoundaryCondition(; transfer_velocity = ScaledTransferVelocity(; schmidt_number = OxygenPolynomialSchmidtNumber()),
                                     water_concentration = ModelConcentration(:O₂),
                                     air_concentration = 9352.7, # mmolO₂/m³
                                     wind_speed = 2) = 
    FluxBoundaryCondition(GasExchange(wind_speed, transfer_velocity, water_concentration, air_concentration, (:O₂, :T), (; T = 1, O₂ = 2));
                          field_dependencies = (:O₂, :T))
end # module