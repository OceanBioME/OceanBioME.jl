include("dependencies_for_runtests.jl")

using Oceananigans, Oceananigans.Units, DataDeps, JLD2, Statistics, Adapt

using Oceananigans.Fields: ConstantField

using OceanBioME: GasExchange, LOBSTER, CarbonChemistry
using OceanBioME.Models.GasExchangeModel: surface_value, CarbonDioxideConcentration, GarciaGordonOxygenSaturation,
                                          CarbonDioxideAirConcentration, PartiallySolubleGas, OxygenSolubility,
                                          MolPerKgPerAtmToMMolPerCubicMPerMicroAtm

using OceanBioME.Models.GasExchangeModel.ScaledGasTransferVelocity: UnitSolubility

using OceanBioME.Models: teos10_polynomial_approximation

using OceanBioME.Models.CarbonChemistryModel: IonicStrength, FF, K0, K1, K2, KB, KW, KS, KF, KP, KSi, KSP_aragonite, KSP_calcite

const year = years = 365days # just for the idealised case below

dd = DataDep(
    "test_data",
    "CODAP-NA (https://essd.copernicus.org/articles/13/2777/2021/) data for testing pCO₂ calculations",
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/CODAP_data.jld2"
)

register(dd)

function test_gas_exchange_model(grid, air_concentration)
    model = NonhydrostaticModel(grid;
                                tracers = (:T, :S),
                                biogeochemistry = LOBSTER(grid; inorganic_carbon = CarbonateSystem()),
                                boundary_conditions = (DIC = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition(; air_concentration)), ))

    set!(model, T = 15.0, S = 35.0, DIC = 2220, Alk = 2500)

    # is everything communicating properly? (can't think of a way to not use allow scalar here)
    value = CUDA.@allowscalar Oceananigans.getbc(model.tracers.DIC.boundary_conditions.top, 1, 1, grid, model.clock, fields(model))

    @test isa(model.tracers.DIC.boundary_conditions.top.condition.func, GasExchange)
    @test ≈(value, -6.8234e-6; atol = 1e-8)
    @test isnothing(time_step!(model, 1.0))

    # multiple carbonate systems

    CO₂_flux1 = 
        CarbonDioxideGasExchangeBoundaryCondition(; 
            water_concentration = CarbonDioxideConcentration(; DIC = :DIC1,
                                                            Alk = :Alk1)
        )
    CO₂_flux2 = 
        CarbonDioxideGasExchangeBoundaryCondition(; 
            water_concentration = CarbonDioxideConcentration(; DIC = :DIC2,
                                                            Alk = :Alk2)
        )


    boundary_conditions = (; DIC1 = FieldBoundaryConditions(top = CO₂_flux1),
                             DIC2 = FieldBoundaryConditions(top = CO₂_flux2))

    model = NonhydrostaticModel(grid;
                                tracers = (:T, :S),
                                biogeochemistry = LOBSTER(grid; inorganic_carbon = CarbonateSystem(2)),
                                boundary_conditions)

    set!(model, T = 15.0, S = 35.0, DIC1 = 2220, Alk1 = 2500, DIC2 = 2221, Alk2 = 2501)

    # maybe this should use BoundaryConditionOperation
    value1 = CUDA.@allowscalar Oceananigans.getbc(model.tracers.DIC1.boundary_conditions.top, 1, 1, grid, model.clock, fields(model))
    value2 = CUDA.@allowscalar Oceananigans.getbc(model.tracers.DIC2.boundary_conditions.top, 1, 1, grid, model.clock, fields(model))

    @test isa(model.tracers.DIC1.boundary_conditions.top.condition.func, GasExchange)
    @test ≈(value1, -6.8234e-6; atol = 1e-8)

    @test isa(model.tracers.DIC2.boundary_conditions.top.condition.func, GasExchange)
    @test ≈(value2, -6.7970e-6; atol = 1e-8)

    @test value1 != value2

    @test isnothing(time_step!(model, 1.0))
    return nothing
end

@testset "pCO₂ values" begin
    # approximatly correct sized pCO₂ from DIC and Alk
    carbon_chem = CarbonChemistry()

    @load datadep"test_data/CODAP_data.jld2" DIC Alk T S pH pCO₂

    pCO₂_results = similar(pCO₂)
    pH_results = similar(pH)

    for (idx, DIC) in enumerate(DIC)
        pCO₂_results[idx] = carbon_chem(; DIC = Float64(DIC), Alk = Float64(Alk[idx]), T = Float64(T[idx]), S = Float64(S[idx]), output = Val(:pCO₂))
        pH_results[idx] = carbon_chem(; DIC, Alk = Alk[idx], T = T[idx], S = S[idx], output = Val(:pHᶠ))
    end

    pCO₂_err = pCO₂ .- pCO₂_results
    pH_err = pH .- pH_results

    # not great not terrible
    @test (mean(pCO₂_err) < 7.1 && std(pCO₂_err) < 9.1) && (mean(pH_err) < 0.01 && std(pH_err) < 0.01)
end

grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 1))

@inline conc_function(x, y, t) = 413.0

conc_field = CenterField(grid)

set!(conc_field, (args...) -> 413)

@testset "Gas exchange coupling" begin
    for air_concentration in [413.1, conc_function, conc_field]
        @info "Testing gas exchange with $(typeof(air_concentration))"
        test_gas_exchange_model(grid, air_concentration)
    end
end

@testset "Carbon chemistry" begin
    carbon_chemistry = CarbonChemistry()

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
    @test ≈(log10(KspA(Tk, S)), -6.1883; atol = 0.0001)
    @test ≈(log10(KspC(Tk, S)), -6.3693; atol = 0.0001)

    @test ≈(KspA.pressure_correction(Tk, P), 1.47866; atol=0.00001)
    @test ≈(KspC.pressure_correction(Tk, P), 1.52962; atol=0.00001)
end

@testset "Gas exchange constants defaults" begin
    for FT in [Float64, Float32]
        CO₂_exchange = CarbonDioxideGasExchangeBoundaryCondition(FT).condition.func
        O₂_exchange = OxygenGasExchangeBoundaryCondition(FT).condition.func

        T = FT(20)

        # values from Wanninkhof, 2014
        @test ≈(CO₂_exchange.transfer_velocity.schmidt_number(T), 668, atol = 1)
        @test ≈(O₂_exchange.transfer_velocity.schmidt_number(T), 568, atol = 1)

        Tk = 25+273.15
        # values from Dickson et. al, 2007
        @test ≈(CO₂_exchange.water_concentration.carbon_chemistry.first_virial_coefficient(Tk), -123.2 * 10^-6, atol=10^-8)
        @test ≈(CO₂_exchange.water_concentration.carbon_chemistry.cross_virial_coefficient(Tk), 22.5 * 10^-6, atol=10^-7)

        T = ConstantField(FT(25))
        S = ConstantField(FT(35))
        DIC = ConstantField(FT(2136.242890518708))
        Alk = ConstantField(FT(2500))
        O₂ = ConstantField(FT(100))

        # value from Dickson et. al, 2007
        pCO₂ = CarbonChemistry(FT)(; DIC = FT(2136.242890518708), Alk = FT(2500), T = FT(25), S = FT(35), output = Val(:pCO₂))
        @test ≈(pCO₂, 350, atol = 0.1)
        @test typeof(pCO₂) == FT

        CO₂ = surface_value(CO₂_exchange.water_concentration, 1, 1, BoxModelGrid(FT), Clock(; time = 0), (; T, S, DIC, Alk))
        @test CO₂ === CarbonChemistry(FT)(; DIC = FT(2136.242890518708), Alk = FT(2500), T = FT(25), S = FT(35), output = Val(:CO₂))
        @test typeof(CO₂) == FT

        pO₂ = surface_value(O₂_exchange.air_concentration, 1, 1, BoxModelGrid(FT), Clock(; time = 0), (; T, S))
        @test ≈(pO₂, 200, atol = 50) # ball park correct
        @test typeof(pO₂) == FT

        # check flux is correct type
        DIC_flux = CO₂_exchange(1, 1, BoxModelGrid(FT), Clock(; time = zero(FT)), (; T, S, DIC, Alk))
        @test typeof(DIC_flux) == FT

        O₂_flux = O₂_exchange(1, 1, BoxModelGrid(FT), Clock(; time = 0), (; T, S, O₂))
        @test typeof(O₂_flux) == FT
    end
end

@testset "Garcia and Gordon (1992) oxygen saturation" begin
    for FT in [Float64, Float32]
        saturation = GarciaGordonOxygenSaturation(FT)

        grid = BoxModelGrid(FT)
        clock = Clock(; time = zero(FT))

        O₂sat(sat, T, S) = surface_value(sat, 1, 1, grid, clock, (T = ConstantField(FT(T)), S = ConstantField(FT(S))))

        # check value from Garcia and Gordon (1992), quoted to six significant figures
        check_value = O₂sat(saturation, 10, 35)

        @test ≈(check_value, 282.015, atol = 5e-4)
        @test typeof(check_value) == FT

        # solubility falls with both temperature and salinity
        @test all(diff([O₂sat(saturation, T, 35) for T in 0:2.5:35]) .< 0)
        @test all(diff([O₂sat(saturation, 10, S) for S in 0:2.5:40]) .< 0)

        # the saturation is exactly linear in the atmospheric pressure, which defaults to 1 atm
        reduced_pressure = GarciaGordonOxygenSaturation(FT; atmospheric_pressure = 0.9)

        @test O₂sat(reduced_pressure, 10, 35) == FT(0.9) * check_value

        # usable as the air concentration of an oxygen gas exchange boundary condition
        exchange = OxygenGasExchangeBoundaryCondition(FT; air_concentration = saturation).condition.func

        O₂ = ConstantField(FT(100))
        T = ConstantField(FT(10))
        S = ConstantField(FT(35))

        flux = exchange(1, 1, grid, clock, (; T, S, O₂))

        @test typeof(flux) == FT
        @test flux < 0 # undersaturated water takes up oxygen

        # GPU compatibility
        @test isbits(saturation)
        @test adapt(Array, saturation) isa GarciaGordonOxygenSaturation
        @test surface_value(adapt(Array, saturation), 1, 1, grid, clock, (; T, S)) == check_value

        # the atmospheric pressure may also be a `Field`
        field_grid = RectilinearGrid(architecture, FT; size = (1, 1, 2), extent = (1, 1, 2))

        pressure_field = CenterField(field_grid)

        set!(pressure_field, 0.9)

        field_pressure = GarciaGordonOxygenSaturation(FT; atmospheric_pressure = pressure_field)

        field_T = CenterField(field_grid)
        field_S = CenterField(field_grid)

        set!(field_T, 10)
        set!(field_S, 35)

        field_value = CUDA.@allowscalar surface_value(field_pressure, 1, 1, field_grid, clock, (T = field_T, S = field_S))

        @test field_value == FT(0.9) * check_value

        adapted = adapt(Array, field_pressure)

        @test CUDA.@allowscalar(surface_value(adapted, 1, 1, field_grid, clock, (T = field_T, S = field_S))) == field_value
    end
end

@testset "Atmospheric pressure on the CO₂ air concentration" begin
    for FT in [Float64, Float32]
        grid = BoxModelGrid(FT)
        clock = Clock(; time = zero(FT))

        T = ConstantField(FT(15))
        S = ConstantField(FT(35))
        DIC = ConstantField(FT(2220))
        Alk = ConstantField(FT(2500))

        model_fields = (; T, S, DIC, Alk)

        # bare (ppmv) air concentrations
        air_concentration = CarbonDioxideAirConcentration(FT)
        reduced_pressure = CarbonDioxideAirConcentration(FT; atmospheric_pressure = 0.9)

        xCO₂ = surface_value(air_concentration, 1, 1, grid, clock, model_fields)

        # the default pressure of 1 atm leaves the mole fraction untouched
        @test xCO₂ === FT(413)

        # the air concentration is exactly linear in the atmospheric pressure
        @test surface_value(reduced_pressure, 1, 1, grid, clock, model_fields) === FT(0.9) * xCO₂

        # only the air side can see the atmospheric pressure: the water-side
        # `CarbonDioxideConcentration` carries none, and rejects one
        carbon_chemistry = CarbonChemistry(FT)

        @test_throws MethodError CarbonDioxideConcentration(FT; carbon_chemistry, air_pressure = FT(0.9))

        # a bare number `air_concentration` still means a mole fraction in ppmv
        default_exchange = CarbonDioxideGasExchangeBoundaryCondition(FT).condition.func
        scalar_exchange = CarbonDioxideGasExchangeBoundaryCondition(FT; air_concentration = 413).condition.func

        @test default_exchange.air_concentration isa CarbonDioxideAirConcentration
        @test !isnothing(default_exchange.air_concentration.solubility) # i.e. it is a concentration
        @test surface_value(default_exchange.air_concentration, 1, 1, grid, clock, model_fields) ===
                surface_value(scalar_exchange.air_concentration, 1, 1, grid, clock, model_fields)
        @test default_exchange(1, 1, grid, clock, model_fields) === scalar_exchange(1, 1, grid, clock, model_fields)

        # the flux is `k (water - air)`, positive out of the ocean, so raising the pressure
        # raises the air term and drives more uptake, and the difference is exactly `k xCO₂ ΔP`
        solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), carbon_chemistry.density_function)

        exchange = CarbonDioxideGasExchangeBoundaryCondition(FT;
                        air_concentration = CarbonDioxideAirConcentration(FT; solubility)).condition.func
        reduced_exchange = CarbonDioxideGasExchangeBoundaryCondition(FT;
                        air_concentration = CarbonDioxideAirConcentration(FT; atmospheric_pressure = 0.9, solubility)).condition.func

        flux = exchange(1, 1, grid, clock, model_fields)
        reduced_flux = reduced_exchange(1, 1, grid, clock, model_fields)

        air_term = surface_value(exchange.air_concentration, 1, 1, grid, clock, model_fields)

        u₁₀ = surface_value(exchange.wind_speed, 1, 1, grid, clock)
        k = exchange.transfer_velocity(u₁₀, FT(15), FT(35))

        @test typeof(flux) == FT
        @test flux < 0             # pCO₂ of 337 μatm under 413 ppmv of air ⇒ uptake
        @test reduced_flux > flux  # less air-side CO₂ ⇒ less uptake
        @test ≈(reduced_flux - flux, k * air_term * FT(0.1); rtol = 100 * eps(FT))

        # GPU compatibility
        @test isbits(air_concentration)
        @test isbits(default_exchange) # the whole assembled exchange, air side, water side and k
        @test adapt(Array, reduced_pressure) isa CarbonDioxideAirConcentration
        @test surface_value(adapt(Array, reduced_pressure), 1, 1, grid, clock, model_fields) ===
                surface_value(reduced_pressure, 1, 1, grid, clock, model_fields)

        # a number, a function and a `Field` pressure all give the same value
        field_grid = RectilinearGrid(architecture, FT; size = (1, 1, 2), extent = (1, 1, 2))

        pressure_field = CenterField(field_grid)

        set!(pressure_field, 0.9)

        field_pressure = CarbonDioxideAirConcentration(FT; atmospheric_pressure = pressure_field)
        function_pressure = CarbonDioxideAirConcentration(FT; atmospheric_pressure = (x, y, t) -> FT(0.9))

        @test CUDA.@allowscalar(surface_value(field_pressure, 1, 1, field_grid, clock, model_fields)) === FT(0.9) * xCO₂
        @test surface_value(function_pressure, 1, 1, field_grid, clock, model_fields) === FT(0.9) * xCO₂

        # `PartiallySolubleGas` adapts a `Field` air concentration
        field_air_concentration = CenterField(field_grid)
        field_T = CenterField(field_grid)
        field_S = CenterField(field_grid)

        set!(field_air_concentration, 9352.7)
        set!(field_T, 10)
        set!(field_S, 35)

        gas = PartiallySolubleGas(FT; air_concentration = field_air_concentration, solubility = OxygenSolubility(FT))

        field_fields = (T = field_T, S = field_S)

        value = CUDA.@allowscalar surface_value(gas, 1, 1, field_grid, clock, field_fields)

        adapted = adapt(Array, gas)

        @test adapted isa PartiallySolubleGas
        @test !(adapted.air_concentration isa Field) # i.e. the `Field` was actually adapted
        @test CUDA.@allowscalar(surface_value(adapted, 1, 1, field_grid, clock, field_fields)) == value
    end
end

# a direct transliteration of MARBL's `ff` (marbl_co2calc_mod.F90:414-423), the authority for the
# `FF` coefficients
function marbl_ff(temp, salt)
    tk = 273.15 + temp
    tk100 = tk * 1e-2
    tk1002 = tk100^2

    arg = -162.8301 + 218.2968 / tk100 + 90.9241 * (log(tk) + log(1e-2)) - 1.47696 * tk1002 +
          salt * (0.025695 - 0.025225 * tk100 + 0.0049867 * tk1002)

    return exp(arg)
end

@testset "CO₂ exchange against MARBL reference values" begin
    for FT in [Float64, Float32]
        grid = BoxModelGrid(FT)
        clock = Clock(; time = zero(FT))

        cc = CarbonChemistry(FT)

        # the reference numbers below are quoted to eight significant figures, which is finer than
        # `Float32` can resolve them
        ref_rtol = max(1e-7, 1000 * eps(FT))

        states = ((T = 10, S = 35, DIC = 2000, Alk = 2000),
                  (T = 15, S = 35, DIC = 2220, Alk = 2500),
                  (T = 25, S = 35, DIC = 2100, Alk = 2350))

        model_fields(s) = (T = ConstantField(FT(s.T)), S = ConstantField(FT(s.S)),
                           DIC = ConstantField(FT(s.DIC)), Alk = ConstantField(FT(s.Alk)))

        # `Val(:CO₂)` is MARBL's `co2star · mass_to_vol`, to machine precision
        # (`co2star = dic H²/(H² + k1 H + k1 k2)`, marbl_co2calc_mod.F90:177)
        for s in states
            args = (; DIC = FT(s.DIC), Alk = FT(s.Alk), T = FT(s.T), S = FT(s.S))

            CO₂ = cc(; args..., output = Val(:CO₂))

            H  = FT(10) ^ -cc(; args..., output = Val(:pHᶠ))
            K₁ = cc.carbonic_acid.K1(FT(s.T) + FT(273.15), FT(s.S))
            K₂ = cc.carbonic_acid.K2(FT(s.T) + FT(273.15), FT(s.S))

            co2star = FT(s.DIC) * H^2 / (H^2 + K₁ * H + K₁ * K₂)

            @test ≈(CO₂, co2star; rtol = 10 * eps(FT))
        end

        # `FF` reproduces the Fortran
        for T in (0, 10, 30), S in (5, 35)
            @test ≈(FF{Float64}()(T + 273.15, Float64(S)), marbl_ff(T, S); rtol = 1e-13)
        end

        # the hand-computed MARBL air term:
        # MARBL's `co2starair · mass_to_vol` = xco2·1e-6 · ff · atmpres · 1e6 ρ_sw, with the
        # constant ρ_sw = 1026 kg/m³ of marbl_constants_mod.F90:59
        marbl_density(args...) = FT(1026)

        marbl_air = CarbonDioxideAirConcentration(FT;
                        solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), marbl_density))

        # the residual is float associativity in the `T² T^2` vs `-1.47696 tk100^2` term, not algebra
        @test ≈(surface_value(marbl_air, 1, 1, grid, clock, model_fields(states[1])),
                413e-6 * marbl_ff(10, 35) * 1.026e6; rtol = max(1e-13, 1000 * eps(FT)))

        @test ≈(surface_value(marbl_air, 1, 1, grid, clock, model_fields(states[1])), 18.295270; rtol = ref_rtol)
        @test ≈(surface_value(marbl_air, 1, 1, grid, clock, model_fields(states[2])), 15.553088; rtol = ref_rtol)

        # our default uses the TEOS-10 density rather than MARBL's constant, which moves the air
        # term by the density ratio and nothing else
        default_exchange = CarbonDioxideGasExchangeBoundaryCondition(FT).condition.func

        for s in states
            ρ = teos10_polynomial_approximation(FT(s.T), FT(s.S))

            @test ≈(surface_value(default_exchange.air_concentration, 1, 1, grid, clock, model_fields(s)),
                    surface_value(marbl_air, 1, 1, grid, clock, model_fields(s)) * ρ / FT(1026); rtol = 1e-6)
        end

        @test ≈(surface_value(default_exchange.air_concentration, 1, 1, grid, clock, model_fields(states[1])),
                18.309978; rtol = ref_rtol)

        # the `K₀`-based air term, for comparison against the converter below
        k0_air_term = K0{FT}()(FT(283.15), FT(35)) * teos10_polynomial_approximation(FT(10), FT(35)) / FT(1e3) * 413

        @test ≈(k0_air_term, 18.608272; rtol = ref_rtol)
        @test k0_air_term > surface_value(default_exchange.air_concentration, 1, 1, grid, clock, model_fields(states[1]))

        # the converter means mmol/m³ per μatm
        converter_air = CarbonDioxideAirConcentration(FT;
            solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(cc.solubility, cc.density_function))

        @test ≈(surface_value(converter_air, 1, 1, grid, clock, model_fields(states[1])),
                k0_air_term; rtol = 10 * eps(FT))

        for s in states
            args = (; DIC = FT(s.DIC), Alk = FT(s.Alk), T = FT(s.T), S = FT(s.S))

            water = surface_value(default_exchange.water_concentration, 1, 1, grid, clock, model_fields(s))

            @test water === cc(; args..., output = Val(:CO₂))
        end

        @test ≈(surface_value(default_exchange.water_concentration, 1, 1, grid, clock, model_fields(states[1])),
                58.937761; rtol = ref_rtol)
        @test ≈(surface_value(default_exchange.water_concentration, 1, 1, grid, clock, model_fields(states[2])),
                12.907474; rtol = ref_rtol)

        # the transfer velocity is a bare piston velocity
        @test default_exchange.transfer_velocity.solubility === UnitSolubility()
        @test default_exchange.transfer_velocity(FT(2), FT(10), FT(35)) ===
                default_exchange.transfer_velocity.base_transfer_velocity(FT(2)) /
                    sqrt(default_exchange.transfer_velocity.schmidt_number(FT(10)) / FT(660))

        # a `carbon_chemistry` is required to build the default `air_concentration`
        @test_throws ArgumentError CarbonDioxideGasExchangeBoundaryCondition(FT; carbon_chemistry = nothing,
                                                                                water_concentration = CarbonDioxideConcentration(FT))
    end
end
