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
       OxygenPolynomialSchmidtNumber,
       GarciaGordonOxygenSaturation,
       CarbonDioxideAirConcentration

using Adapt
using Oceananigans.BoundaryConditions: FluxBoundaryCondition
using Oceananigans.Fields: Center
using Oceananigans.Grids: xnode, ynode

using OceanBioME.Models.CarbonChemistryModel: 
    CarbonChemistry,
    FF,
    silicate_concentration,
    phosphate_concentration

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
                                              transfer_velocity = nothing,
                                              air_concentration = nothing,
                                              wind_speed = 2,
                                              water_concentration = nothing,
                                              silicate_and_phosphate_names = nothing,
                                              kwargs...)

Returns a `FluxBoundaryCondition` for the gas exchange between carbon dioxide dissolved in the water
specified by the `carbon_chemisty` model, and `air_concentration` with `transfer_velocity` (see 
`GasExchangeBoundaryCondition` for details).

By default the exchange is computed on a *concentration* basis, i.e.

    k (‌[CO₂(aq)] - xCO₂ pₐₜₘ ff(T, S) ρ / 10³),

where the water side is the aqueous carbon dioxide concentration in mmol / m³
([`CarbonDioxideConcentration`](@ref) with `output = Val(:CO₂)`), the air side is the dry air
mole fraction converted to a concentration by the Weiss and Price (1980) solubility
([`CarbonDioxideAirConcentration`](@ref)), and `k` is a bare piston velocity.

Passing `water_concentration = CarbonDioxideConcentration(FT; carbon_chemistry, output = Val(:pCO₂))`
instead selects the legacy *partial pressure* basis, in which the solubility is carried by the
`transfer_velocity` and the air concentration is a mole fraction in ppmv. The defaults for
`transfer_velocity` and `air_concentration` are chosen to match whichever basis the
`water_concentration` is on, so the two cannot be mixed by accident.

An `air_concentration` given as a bare number, a function, or a `Field` is interpreted as a dry air
mole fraction in ppmv and is put on the same basis as the `water_concentration`; pass a
[`CarbonDioxideAirConcentration`](@ref) to control that conversion.

`silicate_and_phosphate_names` should either be `nothing`, a `Tuple`` of symbols specifying the name of the silicate
and phosphate tracers, or a `NamedTuple`  of values for the `carbon_chemistry` model.

`kwargs` are passed on to `GasExchangeBoundaryCondition`.

Note: The model always requires `T`, `S`, `DIC`, and `Alk` to be present in the model.
"""
function CarbonDioxideGasExchangeBoundaryCondition(FT = Float64; 
                                                   carbon_chemistry = CarbonChemistry(FT),
                                                   transfer_velocity = nothing,
                                                   air_concentration = nothing,
                                                   wind_speed = 2,
                                                   water_concentration = nothing,
                                                   kwargs...)

    if isnothing(water_concentration)
        water_concentration = CarbonDioxideConcentration(FT; carbon_chemistry)
    elseif !isnothing(carbon_chemistry)
        @warn "Make sure that the `carbon_chemistry` $(carbon_chemistry) is the same as that in `water_concentration` $(water_concentration) (or set it to `nothing`)"
    end

    # these have to be resolved here rather than in the signature since they depend on which basis
    # `water_concentration` is on, which is only known once it has been defaulted
    isnothing(transfer_velocity) && 
        (transfer_velocity = default_carbon_dioxide_transfer_velocity(FT, water_concentration, carbon_chemistry))

    air_concentration = carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry, air_concentration)

    return GasExchangeBoundaryCondition(FT; water_concentration, air_concentration, transfer_velocity, wind_speed, kwargs...)
end

# a bare number, function, or `Field` `air_concentration` has always meant a dry air mole fraction in
# ppmv, so it is wrapped onto whichever basis the `water_concentration` is on rather than being
# compared against a concentration as if it were one
carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry, mole_fraction) = 
    default_carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry; mole_fraction)

carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry, ::Nothing) = 
    default_carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry)

# already carries its own basis
carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry, air_concentration::CarbonDioxideAirConcentration) = 
    air_concentration

const LegacyCarbonDioxideConcentration = CarbonDioxideConcentration{<:Any, <:Any, <:Any, Val{:pCO₂}}

# the concentration basis: MARBL's bare piston velocity, with the solubility on the air side
default_carbon_dioxide_transfer_velocity(FT, water_concentration, carbon_chemistry) = 
    SchmidtScaledTransferVelocity(FT; schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT))

default_carbon_dioxide_air_concentration(FT, water_concentration, carbon_chemistry; mole_fraction = 413) = 
    CarbonDioxideAirConcentration(FT; mole_fraction,
                                      solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), 
                                                                                            carbon_chemistry.density_function))

default_carbon_dioxide_air_concentration(FT, water_concentration, ::Nothing; mole_fraction = 413) = 
    throw(ArgumentError("`CarbonDioxideGasExchangeBoundaryCondition` needs a `carbon_chemistry` to build the " *
                        "default `air_concentration` (it supplies the density), so please pass either a " *
                        "`carbon_chemistry` or an explicit `air_concentration`."))

# the legacy partial pressure basis: the solubility is carried by the transfer velocity instead
default_carbon_dioxide_transfer_velocity(FT, ::LegacyCarbonDioxideConcentration, carbon_chemistry) = 
    SchmidtScaledTransferVelocity(FT; 
        schmidt_number = CarbonDioxidePolynomialSchmidtNumber(FT),
        solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(carbon_chemistry.solubility,
                                                             carbon_chemistry.density_function))

default_carbon_dioxide_transfer_velocity(FT, ::LegacyCarbonDioxideConcentration, ::Nothing) = 
    throw(ArgumentError("`CarbonDioxideGasExchangeBoundaryCondition` needs a `carbon_chemistry` to build the " *
                        "default `transfer_velocity` for the `Val(:pCO₂)` basis (it supplies the solubility " *
                        "and density), so please pass either a `carbon_chemistry` or an explicit " *
                        "`transfer_velocity`."))

default_carbon_dioxide_air_concentration(FT, ::LegacyCarbonDioxideConcentration, carbon_chemistry; mole_fraction = 413) = 
    CarbonDioxideAirConcentration(FT; mole_fraction)

default_carbon_dioxide_air_concentration(FT, ::LegacyCarbonDioxideConcentration, ::Nothing; mole_fraction = 413) = 
    CarbonDioxideAirConcentration(FT; mole_fraction)

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
