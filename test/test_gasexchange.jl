include("dependencies_for_runtests.jl")

using OceanBioME: Boundaries, GasExchange, LOBSTER
using Oceananigans, DataDeps, JLD2, Statistics
using Oceananigans.Units

const year = years = 365days # just for the idealised case below

dd = DataDep(
    "test_data",
    "CODAP-NA (https://essd.copernicus.org/articles/13/2777/2021/) data for testing pCO₂ calculations", 
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/CODAP_data.jld2"
)

register(dd)

function test_gas_exchange_model(grid, air_concentration)
    model = NonhydrostaticModel(; grid, 
                                  tracers = (:T, :S),
                                  biogeochemistry = LOBSTER(; grid, carbonates = true), 
                                  boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; air_concentration, gas = :CO₂)), ))

    @test isa(model.tracers.DIC.boundary_conditions.top.condition.func, GasExchange)

    set!(model, T = 15.0, S = 35.0, DIC = 2220, Alk = 2500)

    # is everything communicating properly? (can't think of a way to not use allow scalar here)
    value = CUDA.@allowscalar Oceananigans.getbc(model.tracers.DIC.boundary_conditions.top, 1, 1, grid, model.clock, fields(model))

    @test value ≈ -0.0003 atol = 0.0001

    # just incase we broke something
    @test isnothing(time_step!(model, 1.0))

    return nothing
end

@testset "Gas exchange values" begin
    # approximatly correct sized pCO₂ from DIC and ALK
    pCO₂_model = Boundaries.OCMIP_default

    @load datadep"test_data/CODAP_data.jld2" DIC Alk T S pH pCO₂
    
    pCO₂_results = similar(pCO₂)

    for (idx, DIC) in enumerate(DIC)
        pCO₂_results[idx] = pCO₂_model(DIC, Alk[idx], T[idx], S[idx])
    end

    pCO₂_err = pCO₂ .- pCO₂_results

    # not great
    @test (mean(pCO₂_err) < 10 && std(pCO₂_err) < 15)
end

grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 1))

@inline conc_function(x, y, t) = 413.0 + 10.0 * sin(t * π / year)

#conc_field = CenterField(grid, indices=(:, :, grid.Nz))
#set!(conc_field, 413.0)
#Oceananigans.fill_halo_regions!(conc_field)

@testset "Gas exchange coupling" begin
    for air_concentration in [413.1, conc_function]#, conc_field]
        @info "Testing with $(typeof(air_concentration))"
        test_gas_exchange_model(grid, air_concentration)
    end
end