using Test
using OceanBioME: Boundaries, GasExchange, LOBSTER
using Oceananigans
using Oceananigans.Grids: xnode, ynode
using Oceananigans.Units

function test_gas_exchange_model(grid, air_concentration)
    PAR = CenterField(grid)
    model = NonhydrostaticModel(;grid, 
                                 tracers = (:T, :S),
                                 biogeochemistry = LOBSTER(;grid, carbonates = true), 
                                 auxiliary_fields=(; PAR),
                                 boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; air_concentration, gas = :CO₂)), ))

    @test isa(model.tracers.DIC.boundary_conditions.top.condition.parameters, GasExchange)

    set!(model, T = 15.0, S = 35.0, DIC = 2220, Alk = 2500)

    # is everything communicating properly?
    @test Oceananigans.getbc(model.tracers.DIC.boundary_conditions.top, 1, 1, grid, model.clock, fields(model)) ≈ -0.0003 atol = 0.0001

    # just incase we broke something
    @test isnothing(time_step!(model, 1.0))

    return nothing
end

@testset "Gas exchange values" begin
    # approximatly correct sized pCO₂ from DIC and ALK
    pCO₂ = Boundaries.pCO₂(2220, 2500, 15 + 273.15, 35.0, 1026.0, 8.0)
    @test pCO₂ ≈ 400 atol = 100 # don't actually care about the value but very easy to mess this up to give wrong order

    # correct sized flux based on pCO₂ = 355 atm and NPP of ∼3mg C/m³/day in top 100m (from Copernicus data) -> ∼0.0003 mmol C/m²/s
    CO₂_exchange_model = GasExchange(;gas = :CO₂)
    @test CO₂_exchange_model.condition.parameters(0.0, 0.0, 0.0, 350.0, 15.0, 35.0) ≈ -0.0003 atol = 0.0001

    # Oxygen exchange (since most gases use a slightly different formulation), must also be same order as NPP*Rd_resp ≈ NPP/0.8 -> ∼0.0004 mmol O₂/m²/s
    O₂_exchange_model = GasExchange(;gas = :O₂)
    @test O₂_exchange_model.condition.parameters(0.0, 0.0, 0.0, 270.0, 15.0, 35.0) ≈ -0.0004 atol = 0.0001
end

@testset "Gas exchange coupling" begin
    grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
    conc_field = CenterField(grid, indices=(:, :, 1))
    conc_field .= 410.0 + 10.0*rand()
    for air_concentration in [410, 410.1, (x, y, t) -> 410.0 + 10.0*sin(t*π/(1year)), conc_field]
        @info "Testing with $(typeof(air_concentration))"
        test_gas_exchange_model(grid, air_concentration)
    end
end