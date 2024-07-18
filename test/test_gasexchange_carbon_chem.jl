include("dependencies_for_runtests.jl")

using OceanBioME: GasExchange, LOBSTER, CarbonChemistry
using Oceananigans, DataDeps, JLD2, Statistics
using Oceananigans.Units

using OceanBioME.Models.CarbonChemistryModel: IonicStrength, K0, K1, K2, KB, KW, KS, KF, KP, KSi

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
    pCO₂_model = CarbonChemistry()

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

@testset "Carbon chemistry" begin
    # vary from default to use all parameters from for testing Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
    # Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

    carbon_chemistry =
        CarbonChemistry(; # Weiss, R.F. (1974, Mar. Chem., 2, 203–215)
                          solubility = K0(-60.2409, 93.4517 * 100, 23.3585, 0.0, 0.023517, -0.023656 / 100, 0.0047036 / 100^2),
                          # Lueke, et. al (2000, Mar. Chem., 70, 105–119;)
                          carbonic_acid = (K1 = K1(61.2172, -3633.86, -9.67770, 0.011555, -0.0001152), 
                                           K2 = K2(-25.9290, -471.78, 0.01781, -0.0001122,  3.16967)),
                          # Perez and Fraga (1987, Mar. Chem., 21, 161–168).
                          fluoride = KF(IonicStrength(),  KS(), -9.68, 874.0, 0.111, 0.0, 0.0))

    S = 35
    Tk = 298.15

    @test ≈(log10(carbon_chemistry.carbonic_acid.K1(Tk, S)), -5.8472; atol=0.0001)
    @test ≈(log10(carbon_chemistry.carbonic_acid.K2(Tk, S)), -8.9660; atol=0.0001)
    @test ≈(log(carbon_chemistry.boric_acid(Tk, S)), -19.7964; atol=0.0001)
    @test ≈(log(carbon_chemistry.water(Tk, S)), -30.434; atol=0.001)
    @test ≈(log(carbon_chemistry.sulfate(Tk, S)), -2.30; atol=0.01)
    @test ≈(log(carbon_chemistry.fluoride(Tk, S)), -6.09; atol=0.01)
    @test ≈(log(carbon_chemistry.phosphoric_acid.KP1(Tk, S)), -3.71; atol=0.01)
    @test ≈(log(carbon_chemistry.phosphoric_acid.KP2(Tk, S)), -13.727; atol=0.001)
    @test ≈(log(carbon_chemistry.phosphoric_acid.KP3(Tk, S)), -20.24; atol=0.01)
    @test ≈(log(carbon_chemistry.silicic_acid(Tk, S)), -21.61; atol=0.01)
    @test ≈(log(carbon_chemistry.solubility(Tk, S)), -3.5617; atol=0.0001)
end