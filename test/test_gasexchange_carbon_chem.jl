include("dependencies_for_runtests.jl")

using OceanBioME: GasExchange, LOBSTER, CarbonChemistry
using Oceananigans, DataDeps, JLD2, Statistics
using Oceananigans.Units

using OceanBioME.Models.CarbonChemistryModel: IonicStrength, K0, K1, K2, KB, KW, KS, KF, KP, KSi, KSP_aragonite, KSP_calcite

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
                                  boundary_conditions = (DIC = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition(; air_concentration)), ))

    set!(model, T = 15.0, S = 35.0, DIC = 2220, Alk = 2500)

    # is everything communicating properly? (can't think of a way to not use allow scalar here)
    value = CUDA.@allowscalar Oceananigans.getbc(model.tracers.DIC.boundary_conditions.top, 1, 1, grid, model.clock, fields(model))

    return isa(model.tracers.DIC.boundary_conditions.top.condition.func, GasExchange)&&≈(value, -0.0002; atol = 0.0001)&&isnothing(time_step!(model, 1.0))
end

@testset "pCO₂ values" begin
    # approximatly correct sized pCO₂ from DIC and Alk
    fCO₂_model = CarbonChemistry()

    @load datadep"test_data/CODAP_data.jld2" DIC Alk T S pH pCO₂
    
    fCO₂_results = similar(pCO₂)
    pH_results = similar(pH)

    for (idx, DIC) in enumerate(DIC)
        fCO₂_results[idx] = fCO₂_model(; DIC, Alk = Alk[idx], T = T[idx], S = S[idx])
        pH_results[idx] = fCO₂_model(; DIC, Alk = Alk[idx], T = T[idx], S = S[idx], return_pH = true)
    end

    fCO₂_err = pCO₂ .- fCO₂_results
    pH_err = pH .- pH_results

    # not great not terrible
    @test (mean(fCO₂_err) < 9 && std(fCO₂_err) < 10) && (mean(pH_err) < 0.01 && std(pH_err) < 0.01)
end

grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 1))

@inline conc_function(x, y, t) = 413.0 + 10.0 * sin(t * π / year)

@testset "Gas exchange coupling" begin
    for air_concentration in [413.1, conc_function]
        @info "Testing with $(typeof(air_concentration))"
        if air_concentration == conc_function && isa(architecture, GPU) # not sure why, adapt won't work for it
            @test_broken test_gas_exchange_model(grid, air_concentration)
        else
            @test test_gas_exchange_model(grid, air_concentration)
        end
    end
end

@testset "Carbon chemistry" begin
    # vary from default to use all parameters from for testing Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
    # Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

    carbon_chemistry =
        CarbonChemistry(; # Weiss, R.F. (1974, Mar. Chem., 2, 203–215)
                          solubility = K0(-60.2409, 93.4517 * 100, 23.3585, 0.0, 0.023517, -0.023656 / 100, 0.0047036 / 100^2),
                          # Lueke, et. al (2000, Mar. Chem., 70, 105–119;)
                          carbonic_acid = (K1 = K1(constant=61.2172, inverse_T=-3633.86, log_T=-9.67770, S=0.011555, S²=-0.0001152), 
                                           K2 = K2(constant=-25.9290, inverse_T=-471.78, log_T=3.16967, S=0.01781, S²=-0.0001122)),
                          # Perez and Fraga (1987, Mar. Chem., 21, 161–168).
                          fluoride = KF(constant=-9.68, inverse_T=874.0, sqrt_S=0.111, log_S=0.0, log_S_KS=0.0))

    # test conditions
    S = 35
    Tk = 298.15
    P = 300
    
    # values from Dickson et. al, 2007
    @test ≈(log(carbon_chemistry.solubility(Tk, S)), -3.5617; atol=0.0001)
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

    # values from Zeebe & Wolf-Gladrow, 2001
    @test ≈(carbon_chemistry.carbonic_acid.K1.pressure_correction(Tk, P), 1.30804; atol=0.00001)
    @test ≈(carbon_chemistry.carbonic_acid.K2.pressure_correction(Tk, P), 1.21341; atol=0.00001)
    @test ≈(carbon_chemistry.boric_acid.pressure_correction(Tk, P), 1.38024; atol=0.00001)
    @test ≈(carbon_chemistry.water.pressure_correction(Tk, P), 1.23784; atol=0.00001)
    @test ≈(carbon_chemistry.sulfate.pressure_correction(Tk, P), 1.21844; atol=0.00001)
    @test ≈(carbon_chemistry.fluoride.pressure_correction(Tk, P), 1.13151; atol=0.00001)
    @test ≈(carbon_chemistry.phosphoric_acid.KP1.pressure_correction(Tk, P), 1.14852; atol=0.00001)
    @test ≈(carbon_chemistry.phosphoric_acid.KP2.pressure_correction(Tk, P), 1.27298; atol=0.00001)
    @test ≈(carbon_chemistry.phosphoric_acid.KP3.pressure_correction(Tk, P), 1.32217; atol=0.00001)

    # Calcite and aragonite solubility
    KspA = KSP_aragonite()
    KspC = KSP_calcite()

    # Zeebe & Wolf-Gladrow, 2001, Appendix A
    @test ≈(log(KspA(Tk, S)), -6.1883; atol = 0.0001)
    @test ≈(log(KspC(Tk, S)), -6.3693; atol = 0.0001)

    @test ≈(KspA.pressure_correction(Tk, P), 1.47866; atol=0.00001)
    @test ≈(KspC.pressure_correction(Tk, P), 1.52962; atol=0.00001)
end

@testset "Gas exchange constants defaults" begin
    CO₂_exchange = CarbonDioxideGasExchangeBoundaryCondition().condition.func
    O₂_exchange = OxygenGasExchangeBoundaryCondition().condition.func

    T = 20

    # values from Wanninkhof, 2014
    @test ≈(CO₂_exchange.transfer_velocity.schmidt_number(T), 668, atol = 1)
    @test ≈(O₂_exchange.transfer_velocity.schmidt_number(T), 568, atol = 1)

    Tk = 25+273.15
    # values from Dickson et. al, 2007
    @test ≈(CO₂_exchange.water_concentration.first_virial_coefficient(Tk), -123.2 * 10^-6, atol=10^-8)
    @test ≈(CO₂_exchange.water_concentration.cross_virial_coefficient(Tk), 22.5 * 10^-6, atol=10^-7)

    T = 25
    S = 35
    p = 1
    DIC = 2136.242890518708
    Alk = 2500
    # value from Dickson et. al, 2007
    @test ≈(CO₂_exchange.water_concentration(0, 0, 0, T, S, DIC, Alk), 350, atol = 0.1)
end