include("dependencies_for_runtests.jl")

using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: ConstantField, FunctionField

using OceanBioME.Models.PISCESModel: SimpleIron, NitrateAmmonia

const PISCES_INITIAL_VALUES = (P = 0.5, PChl = 0.02, PFe = 0.005,
                               D = 0.1, DChl = 0.004, DFe = 0.001, DSi = 0.01,
                               Z = 0.1, M = 0.7,
                               DOC = 2.1, 
                               POC = 7.8, SFe = 0.206, 
                               GOC = 38, BFe = 1.1, PSi = 0.1, CaCO₃ = 10^-10,
                               NO₃ = 2.3, NH₄ = 0.9, PO₄ = 0.6, Fe = 0.13, Si = 8.5,
                               DIC = 2205, Alk = 2566, O₂ = 317,
                               T = 10, S = 35)

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

function test_PISCES_conservation() # only on CPU please
    @info "Testing PISCES element conservation (C, Fe, P, Si, N)"

    validation_warning = "This implementation of PISCES is in early development and has not yet been validated against the operational version"

    grid = BoxModelGrid(; z = -5)

    PAR₁ = ConstantField(100)
    PAR₂ = ConstantField(100)
    PAR₃ = ConstantField(100)
    PAR  = ConstantField(300)

    mixed_layer_depth = ConstantField(-10)
    euphotic_depth  = ConstantField(-10)
    euphotic_depth  = ConstantField(-10)
    mean_mixed_layer_vertical_diffusivity = ConstantField(1)
    mean_mixed_layer_light = ConstantField(300)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = @test_warn validation_warning PISCES(; grid, 
                                                             sinking_speeds = (POC = 0, GOC = 0), 
                                                             light_attenuation, 
                                                             mixed_layer_depth, 
                                                             euphotic_depth,
                                                             mean_mixed_layer_light,
                                                             mean_mixed_layer_vertical_diffusivity,
                                                             # turn off permanent iron removal and nitrogen fixaiton
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

    # double precision floats are only valid to 17 bits so this tollerance is actually good
    @test isapprox(total_carbon_tendencies, 0, atol = 10^-20) 
    @test isapprox(total_iron_tendencies, 0, atol = 10^-21) 
    @test isapprox(total_silicon_tendencies, 0, atol = 10^-30) 
    @test isapprox(total_phosphate_tendencies, 0, atol = 10^-22) 
    @test isapprox(total_nitrogen_tendencies, 0, atol = 10^-21) 

    return nothing
end

@inline total_light(z) = 3light(z)
@inline light(z) = ifelse(z <= 0, exp(z/10), 2-exp(-z/10)) # so we get a value boundary condition like normal PAR fields
@inline κ_func(z) = ifelse(z > -25, 2, 1)

function test_PISCES_update_state(arch)
    @info "Testing PISCES auxiliary field computation and timestepping"
    # TODO: implement and test mixed layer depth computaiton elsewhere

    grid = RectilinearGrid(arch, topology = (Flat, Flat, Bounded), size = (10, ), extent = (100, ))

    PAR₁ = PAR₂ = PAR₃ = FunctionField{Center, Center, Center}(light, grid)
    PAR  = FunctionField{Center, Center, Center}(total_light, grid)

    mixed_layer_depth = ConstantField(-25)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = PISCES(; grid, 
                               light_attenuation,
                               mixed_layer_depth)

    closure = ScalarDiffusivity(ν = 1e-2, κ = FunctionField{Center, Center, Center}((z) -> ifelse(z > -25, 2, 1), grid))

    model = NonhydrostaticModel(; grid, biogeochemistry, closure) # this updates the biogeochemical state

    # checked and at very high resolution this converges to higher tollerance
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

    mixed_layer_depth = ConstantField(-25)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = PISCES(; grid, 
                               light_attenuation,
                               mixed_layer_depth,
                               scale_negatives = true)

    model = NonhydrostaticModel(; grid, biogeochemistry)

    set!(model, P = -1, D = 1, Z = 1, M = 1, DOC = 1, POC = 1, GOC = 1, DIC = 1, CaCO₃ = 1, PO₄ = 1)

    # got rid of the negative
    @test on_architecture(CPU(), interior(model.tracers.P, 1, 1, 1))[1] == 0

    # correctly conserved mass
    @test all(map(t -> on_architecture(CPU(), interior(t, 1, 1, 1))[1] ≈ 7/8, model.tracers[(:D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃)]))

    # didn't touch the others
    @test on_architecture(CPU(), interior(model.tracers.PO₄, 1, 1, 1))[1] == 1

    # failed to scale silcate since nothing else in its group was available
    set!(model, Si = -1, DSi = 0.1)
    
    @test isnan(on_architecture(CPU(), interior(model.tracers.DSi, 1, 1, 1))[1])

    # this is actually going to cause failures in conserving other groups because now we'll have less carbon etc from the M and Z...
    set!(model, Fe = -1, Z = 1000, M = 0)

    @test on_architecture(CPU(), interior(model.tracers.Fe, 1, 1, 1))[1] == 0
    @test on_architecture(CPU(), interior(model.tracers.Z, 1, 1, 1))[1] ≈ 900
end
#=
@testset "PISCES" begin
    if architecture isa CPU
        test_PISCES_conservation()
        #test_PISCES_box_model() #TODO
    end

    test_PISCES_update_state(architecture)

    test_PISCES_negativity_protection(architecture)

    #test_PISCES_setup(grid) # maybe should test everything works with all the different bits???
end=#