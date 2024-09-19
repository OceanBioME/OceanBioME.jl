include("dependencies_for_runtests.jl")

using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: ConstantField, FunctionField

using OceanBioME.Models.PISCESModel: SimpleIron, NitrateAmmonia

const PISCES_INITIAL_VALUES = (P = 7, D = 7, Z = 0.5, M = 0.5, PChl = 1.5, DChl = 1.5, PFe = 35e-3, DFe = 35e-3, DSi = 1,
                               NO₃ = 6, NH₄ = 1, PO₄ = 1, Fe = 1, Si = 7, CaCO₃ = 10^-3,
                               DIC = 2200, Alk = 2400, O₂ = 240, T = 10, S = 35)


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
    validation_warning = "This implementation of PISCES is in early development and has not yet been validated"

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
                                                             iron = SimpleIron(0),
                                                             nitrogen = NitrateAmmonia(maximum_fixation_rate = 0.0))

    model = BoxModel(; grid, biogeochemistry)

    # checks model works with zero values
    time_step!(model, 1.0)

    # and that they all return zero
    @test all([all(!(name in (:T, :S)) | (Array(interior(values))[1] .== 0)) for (name, values) in pairs(model.fields)]) 

    # set some semi-realistic conditions and check conservation
    set_PISCES_initial_values!(model.fields)

    time_step!(model, 1.0)

    conserved_tracers = OceanBioME.conserved_tracers(biogeochemistry)

    total_carbon_tendencies = sum(map(f -> value(f), model.timestepper.Gⁿ[conserved_tracers.carbon]))
    total_iron_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.iron.tracers, conserved_tracers.iron.scalefactors)])
    total_silicon_tendencies = sum(map(f -> value(f), model.timestepper.Gⁿ[conserved_tracers.silicon]))
    total_phosphate_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.phosphate.tracers, conserved_tracers.phosphate.scalefactors)])
    total_nitrogen_tendencies = sum([value(model.timestepper.Gⁿ[n]) * sf for (n, sf) in zip(conserved_tracers.nitrogen.tracers, conserved_tracers.nitrogen.scalefactors)])

    # should these be exactly zero?
    @test isapprox(total_carbon_tendencies, 0, atol = 10^-20) 
    @test isapprox(total_iron_tendencies, 0, atol = 10^-21) 
    @test isapprox(total_silicon_tendencies, 0, atol = 10^-21) 
    @test isapprox(total_phosphate_tendencies, 0, atol = 10^-21) 
    @test isapprox(total_nitrogen_tendencies, 0, atol = 10^-21) 

    return nothing
end

total_light(z) = 3light(z)
light(z) = ifelse(z <= 0, exp(z/10), 2-exp(-z/10)) # so we get a value boundary condition like normal PAR fields

function test_PISCES_update_state(arch)
    # TODO: implement and test mixed layer depth computaiton elsewhere

    grid = RectilinearGrid(arch, topology = (Flat, Flat, Bounded), size = (10, ), extent = (100, ))

    PAR₁ = PAR₂ = PAR₃ = FunctionField{Center, Center, Center}(light, grid)
    PAR  = FunctionField{Center, Center, Center}(total_light, grid)

    mixed_layer_depth = ConstantField(-25)

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

    biogeochemistry = PISCES(; grid, 
                               light_attenuation,
                               mixed_layer_depth)

    closure = ScalarDiffusivity(ν = nothing, κ = FunctionField{Center, Center, Center}((z) -> ifelse(z > -25, 2, 1), grid))

    model = NonhydrostaticModel(; grid, biogeochemistry, closure) # this updates the biogeochemical state

    # checked and at very high resolution this converges to higher tollerance
    @test isapprox(biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_light[1, 1, 1], 3 * 10 / 25 * (1 - exp(-25/10)), rtol = 0.1)
    @test isapprox(biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_vertical_diffusivity[1, 1, 1], 2, rtol = 0.1)

    # test should be elsewhere
    @test biogeochemistry.underlying_biogeochemistry.euphotic_depth[1, 1, 1] ≈ -10 * log(1000)
end

@testset "PISCES" begin
    if architecture isa CPU
        @info "Testing PISCES element conservation (C, Fe, P, Si)"
        test_PISCES_conservation()
    end

    @info "Testing PISCES auxiliary field computation"
    test_PISCES_update_state(architecture)

    
    #test_PISCES_setup(grid)
end