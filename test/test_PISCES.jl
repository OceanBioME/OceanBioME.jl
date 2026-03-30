include("dependencies_for_runtests.jl")
architecture = CPU()

using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: ConstantField, FunctionField

using OceanBioME.Models.PISCESModel: PISCES, SimpleIron, NitrateAmmonia
using OceanBioME.Models.PISCESModel.Iron: IronTendencyArgs, iron_tendency, ligand_aggregation, ligand_concentration, free_iron

const PISCES_INITIAL_VALUES = (P = 0.5, PChl = 0.02, PFe = 0.005,
                               D = 0.1, DChl = 0.004, DFe = 0.001, DSi = 0.01,
                               Z = 0.1, M = 0.7,
                               DOC = 2.1,
                               POC = 7.8, SFe = 0.206,
                               GOC = 38, BFe = 1.1, PSi = 0.1, CaCO₃ = 10^-10,
                               NO₃ = 2.3, NH₄ = 0.9, PO₄ = 0.6, Fe = 0.13, Si = 8.5,
                               DIC = 2205.0, Alk = 2566.0, O₂ = 317.0,
                               T = 10.0, S = 35.0)

function set_PISCES_initial_values!(tracers; values = PISCES_INITIAL_VALUES)
    for (name, field) in pairs(tracers)
        if name in keys(values)
            set!(field, values[name])
        else
            set!(field, 0)
        end
    end

    return nothing
end

value(field; indices = (1, 1, 1)) = on_architecture(CPU(), interior(field, indices...))[1]

const PISCES_SCALAR_IRON_TEST_VALUES = merge(PISCES_INITIAL_VALUES, (;
    O₂ = 1.0,
    T = 10.0,
    z = -5.0,
    zₘₓₗ = -10.0,
    zₑᵤ = -10.0,
    Si′ = 7.5,
    background_shear = 0.01,
    mixed_layer_shear = 1.0,
    sinking_flux = 0.0,
    sinking_iron_flux = 0.0,
))

all_scalar_numbers(result::Number) = true
all_scalar_numbers(result::Tuple) = all(all_scalar_numbers, result)
all_scalar_numbers(result) = false

function scalar_iron_tendency(bgc; values = PISCES_SCALAR_IRON_TEST_VALUES)
    pom = bgc.particulate_organic_matter
    (aggregation_parameter_1,
     aggregation_parameter_2,
     aggregation_parameter_3,
     aggregation_parameter_4,
     aggregation_parameter_5) = bgc.dissolved_organic_matter.aggregation_parameters

    nano = bgc.phytoplankton.nano
    diatoms = bgc.phytoplankton.diatoms
    zoo = bgc.zooplankton
    micro = zoo.micro
    meso = zoo.meso

    inputs = IronTendencyArgs(
        values.Fe,
        values.DOC,
        values.T,
        values.POC,
        values.GOC,
        values.SFe,
        values.CaCO₃,
        values.PSi,
        values.P,
        values.PChl,
        values.PFe,
        values.D,
        values.DChl,
        values.DFe,
        values.Z,
        values.M,
        values.NH₄,
        values.NO₃,
        values.PO₄,
        values.Si,
        values.O₂,
        values.z,
        values.zₘₓₗ,
        values.zₑᵤ,
        values.Si′,
        values.background_shear,
        values.mixed_layer_shear,
        values.sinking_flux,
        values.sinking_iron_flux,
        bgc.first_anoxia_threshold,
        bgc.second_anoxia_threshold,
        pom.minimum_iron_scavenging_rate,
        pom.load_specific_iron_scavenging_rate,
        pom.base_breakdown_rate,
        pom.temperature_sensitivity,
        pom.maximum_iron_ratio_in_bacteria,
        pom.iron_half_saturation_for_bacteria,
        pom.bacterial_iron_uptake_efficiency,
        pom.maximum_bacterial_growth_rate,
        bgc.dissolved_organic_matter.aggregation_parameters[1],
        bgc.dissolved_organic_matter.aggregation_parameters[2],
        bgc.dissolved_organic_matter.aggregation_parameters[3],
        bgc.dissolved_organic_matter.aggregation_parameters[4],
        bgc.dissolved_organic_matter.aggregation_parameters[5],
        zoo.microzooplankton_bacteria_concentration,
        zoo.mesozooplankton_bacteria_concentration,
        zoo.maximum_bacteria_concentration,
        zoo.bacteria_concentration_depth_exponent,
        zoo.doc_half_saturation_for_bacterial_activity,
        zoo.nitrate_half_saturation_for_bacterial_activity,
        zoo.ammonia_half_saturation_for_bacterial_activity,
        zoo.phosphate_half_saturation_for_bacterial_activity,
        zoo.iron_half_saturation_for_bacterial_activity,
        nano.exudated_fraction,
        nano.maximum_iron_ratio,
        nano.half_saturation_for_iron_uptake,
        nano.threshold_for_size_dependency,
        nano.size_ratio,
        nano.nutrient_limitation.minimum_ammonium_half_saturation,
        nano.nutrient_limitation.minimum_nitrate_half_saturation,
        nano.nutrient_limitation.optimal_iron_quota,
        nano.growth_rate.base_growth_rate,
        nano.growth_rate.temperature_sensitivity,
        diatoms.exudated_fraction,
        diatoms.maximum_iron_ratio,
        diatoms.half_saturation_for_iron_uptake,
        diatoms.threshold_for_size_dependency,
        diatoms.size_ratio,
        diatoms.nutrient_limitation.minimum_ammonium_half_saturation,
        diatoms.nutrient_limitation.minimum_nitrate_half_saturation,
        diatoms.nutrient_limitation.optimal_iron_quota,
        diatoms.growth_rate.base_growth_rate,
        diatoms.growth_rate.temperature_sensitivity,
        micro.iron_ratio,
        micro.non_assimilated_fraction,
        micro.maximum_grazing_rate,
        micro.temperature_sensitivity,
        micro.food_preferences.P,
        micro.food_preferences.D,
        micro.food_preferences.Z,
        micro.food_preferences.POC,
        micro.specific_food_threshold_concentration,
        micro.grazing_half_saturation,
        micro.food_threshold_concentration,
        micro.minimum_growth_efficiency,
        micro.maximum_flux_feeding_rate,
        meso.iron_ratio,
        meso.non_assimilated_fraction,
        meso.maximum_grazing_rate,
        meso.temperature_sensitivity,
        meso.food_preferences.P,
        meso.food_preferences.D,
        meso.food_preferences.Z,
        meso.food_preferences.POC,
        meso.specific_food_threshold_concentration,
        meso.grazing_half_saturation,
        meso.food_threshold_concentration,
        meso.minimum_growth_efficiency,
        meso.maximum_flux_feeding_rate,
        meso.quadratic_mortality,
    )

    return iron_tendency(bgc.iron, inputs)
end

scalar_iron_test_cases(bgc; values = PISCES_SCALAR_IRON_TEST_VALUES) = (
    (; name = "iron_tendency",
       values,
       call = () -> scalar_iron_tendency(bgc; values)),
)

function test_PISCES_scalar_iron_functions()
    @info "Testing scalar PISCES iron functions"

    validation_warning = "This implementation of PISCES is in early development and has not yet been validated against the operational version"

    grid = BoxModelGrid(; z = -5)
    bgc = (@test_warn validation_warning PISCES(; grid)).underlying_biogeochemistry

    for case in scalar_iron_test_cases(bgc)
        @testset "$(case.name)" begin
            result = case.call()

            @test all_scalar_numbers(result)
        end
    end

    return nothing
end

function test_PISCES_conservation() # only on CPU please
    @info "Testing PISCES element conservation (C, Fe, P, Si, N)"

    validation_warning = "This implementation of PISCES is in early development and has not yet been validated against the operational version"

    grid = BoxModelGrid(; z = -5)

    PAR₁ = ConstantField(100.0)
    PAR₂ = ConstantField(100.0)
    PAR₃ = ConstantField(100.0)
    PAR  = ConstantField(300.0)

    mixed_layer_depth = ConstantField(-10.0)
    euphotic_depth  = ConstantField(-10.0)
    mean_mixed_layer_vertical_diffusivity = ConstantField(1.0)
    mean_mixed_layer_light = ConstantField(300.0)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = @test_warn validation_warning PISCES(; grid,
                                                             sinking_speeds = (POC = 0, GOC = 0),
                                                             light_attenuation,
                                                             mixed_layer_depth,
                                                             euphotic_depth,
                                                             mean_mixed_layer_light,
                                                             mean_mixed_layer_vertical_diffusivity,
                                                             # turn off permanent iron removal and nitrogen fixation
                                                             iron = SimpleIron(excess_scavenging_enhancement = 0.0),
                                                             nitrogen = NitrateAmmonia(maximum_fixation_rate = 0.0))

    model = BoxModel(; grid, biogeochemistry)

    # checks model works with zero values
    time_step!(model, 1.0)

    # and that they all return zero
    @test all([all(!(name in (:T, :S)) | (Array(interior(values))[1] .== 0)) for (name, values) in pairs(model.fields)])

    # set some semi-realistic conditions and check conservation
    set_PISCES_initial_values!(model.fields)

    time_step!(model, 1.0)

    conserved_tracers = OceanBioME.conserved_tracers(biogeochemistry; ntuple = true)

    total_carbon_tendencies = sum(map(f -> value(f), model.timestepper.Gⁿ[conserved_tracers.carbon]))
    total_iron_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.iron.tracers, conserved_tracers.iron.scalefactors)])
    total_silicon_tendencies = sum(map(f -> value(f), model.timestepper.Gⁿ[conserved_tracers.silicon]))
    total_phosphate_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.phosphate.tracers, conserved_tracers.phosphate.scalefactors)])
    total_nitrogen_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.nitrogen.tracers, conserved_tracers.nitrogen.scalefactors)])

    # double precision floats are only valid to 17 bits so this tolerance is actually good
    @test isapprox(total_carbon_tendencies, 0, atol = 10^-20)
    @test isapprox(total_iron_tendencies, 0, atol = 10^-21)
    @test isapprox(total_silicon_tendencies, 0, atol = 10^-21)
    @test isapprox(total_phosphate_tendencies, 0, atol = 10^-22)
    @test isapprox(total_nitrogen_tendencies, 0, atol = 10^-21)

    return nothing
end

@inline total_light(z) = 3light(z)
@inline light(z) = ifelse(z <= 0, exp(z/10), 2-exp(-z/10)) # so we get a value boundary condition like normal PAR fields
@inline κ_func(z) = ifelse(z > -25, 2, 1)

function test_PISCES_update_state(arch)
    @info "Testing PISCES auxiliary field computation and time stepping"
    # TODO: implement and test mixed layer depth computation elsewhere

    grid = RectilinearGrid(arch, topology = (Flat, Flat, Bounded), size = (10, ), extent = (100, ))

    PAR₁ = PAR₂ = PAR₃ = FunctionField{Center, Center, Center}(light, grid)
    PAR  = FunctionField{Center, Center, Center}(total_light, grid)

    mixed_layer_depth = ConstantField(-25.0)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = PISCES(; grid,
                               light_attenuation,
                               mixed_layer_depth)

    closure = ScalarDiffusivity(ν = 1e-2, κ = FunctionField{Center, Center, Center}((z) -> ifelse(z > -25, 2, 1), grid))

    model = NonhydrostaticModel(grid; biogeochemistry, closure) # this updates the biogeochemical state

    # checked and at very high resolution this converges to higher tolerance
    @test isapprox(on_architecture(CPU(), biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_light)[1, 1, 1], 3 * 10 / 25 * (1 - exp(-25/10)), rtol = 0.1)
    @test isapprox(on_architecture(CPU(), biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_vertical_diffusivity)[1, 1, 1], 2, rtol = 0.1)

    # test should be elsewhere
    @test on_architecture(CPU(), biogeochemistry.underlying_biogeochemistry.euphotic_depth)[1, 1, 1] ≈ -10 * log(1000)

    @test_nowarn time_step!(model, 1)
end

function test_PISCES_negativity_protection(arch)
    @info "Testing PISCES negativity protection"
    grid = RectilinearGrid(arch, topology = (Flat, Flat, Bounded), size = (10, ), extent = (100, ))

    PAR₁ = PAR₂ = PAR₃ = FunctionField{Center, Center, Center}(light, grid)
    PAR  = FunctionField{Center, Center, Center}(total_light, grid)

    mixed_layer_depth = ConstantField(-25.0)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = PISCES(; grid,
                               light_attenuation,
                               mixed_layer_depth,
                               scale_negatives = true)

    model = NonhydrostaticModel(grid; biogeochemistry)

    set!(model, P = -1, D = 1, Z = 1, M = 1, DOC = 1, POC = 1, GOC = 1, DIC = 1, CaCO₃ = 1, PO₄ = 1)

    # got rid of the negative
    @test on_architecture(CPU(), interior(model.tracers.P, 1, 1, 1))[1] == 0

    # correctly conserved mass
    @test all(map(t -> on_architecture(CPU(), interior(t, 1, 1, 1))[1] ≈ 7/8, model.tracers[(:D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃)]))

    # didn't touch the others
    @test on_architecture(CPU(), interior(model.tracers.PO₄, 1, 1, 1))[1] == 1

    # failed to scale silicate since nothing else in its group was available
    set!(model, Si = -1, DSi = 0.1)

    @test isnan(on_architecture(CPU(), interior(model.tracers.DSi, 1, 1, 1))[1])

    # this is actually going to cause failures in conserving other groups because now we'll have less carbon etc from the M and Z...
    set!(model, Fe = -1, Z = 1000, M = 0)

    @test on_architecture(CPU(), interior(model.tracers.Fe, 1, 1, 1))[1] == 0
    @test on_architecture(CPU(), interior(model.tracers.Z, 1, 1, 1))[1] ≈ 900
end

@testset "PISCES" begin
    if architecture isa CPU
        test_PISCES_conservation()
        #test_PISCES_box_model() #TODO
    end

    test_PISCES_update_state(architecture)

    test_PISCES_negativity_protection(architecture)

    if architecture isa CPU
        test_PISCES_scalar_iron_functions()
    end

    #test_PISCES_setup(grid) # maybe should test everything works with all the different bits???
end
