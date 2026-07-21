include("dependencies_for_runtests.jl")

using OceanBioME: conserved_tracers
using OceanBioME.Models.NutrientsPlanktonDetritusModels: InstantRemineralisationDetritus, 
                                                         CarbonNitrogenDissolvedParticulate
using Oceananigans, CUDA, Random

Random.seed!(42)

grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

test_models = []

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

using CUDA
using OceanBioME: conserved_tracers

function get_conservation_values(bgc, fields)
    tracer_groups = conserved_tracers(bgc)

    conserved_values = NamedTuple()

    for (group_name, group) in pairs(tracer_groups)
        value = 0.0

        for (name, sf) in pairs(group)
            value += sf * (CUDA.@allowscalar first(Field(Integral(fields[name]))))
        end

        conserved_values = merge(conserved_values, NamedTuple{(group_name,)}((value, )))
    end

    return conserved_values
end

# everything should be fine as long as everythings positive right???
set_default!(model) = [set!(tracer, rand()*10) for tracer in model.tracers]

function check_stepped(model)
    set_default!(model)
    initial_values = NamedTuple(name => CUDA.@allowscalar model.tracers[name][1, 1, 1]
                                for name in keys(model.tracers))
    time_step!(model, 1)
    for (name, initial) in pairs(initial_values)
        if !iszero(initial)
            @test CUDA.@allowscalar model.tracers[name][1, 1, 1] != initial
        end
    end

    return nothing
end

function check_conservations(model, n_timesteps = 100)
    set_default!(model)

    initial_values = get_conservation_values(model.biogeochemistry, model.tracers)

    for _ in 1:n_timesteps
        time_step!(model, 1)
    end

    final_values = get_conservation_values(model.biogeochemistry, model.tracers)

    for (name, value) in pairs(initial_values)
        @info "Checking conservation of $(name)"
        @test isapprox(final_values[name], value, atol = n_timesteps*eps(value))
    end
end

# need to do at least one test of detritus and dissolved particulate with advection on
nutrients_options = (Nutrients(; nitrogen = OceanBioME.N),
                     Nutrients(; nitrogen = NitrateAmmonia()),
                     Nutrients(; phosphate = OceanBioME.PO₄),
                     Nutrients(; iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = OceanBioME.N, phosphate = OceanBioME.PO₄),
                     Nutrients(; nitrogen = OceanBioME.N, iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = NitrateAmmonia(), phosphate = OceanBioME.PO₄),
                     Nutrients(; nitrogen = NitrateAmmonia(), iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = OceanBioME.N, phosphate = OceanBioME.PO₄, iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = NitrateAmmonia(), phosphate = OceanBioME.PO₄, iron = OceanBioME.Fe))

detritus_options = (InstantRemineralisationDetritus(),
                    Detritus(grid), 
                    DissolvedParticulate(grid, :DOP, :POP), 
                    DissolvedParticulate(grid),
                    CarbonNitrogenDissolvedParticulate(grid))

# TODO: test multiple separately
inorganic_carbon_options = (nothing, 
                            CarbonateSystem(),)#, CarbonateSystem(2))

oxygen_options = (nothing, 
                  Oxygen(),)

@testset "Elemental conservations" begin 
    # maybe this is pointless but would like to keep Abiotic working for testing
    # (abiotic is going to result in everything being remineralised since theres no path to make detritus)
    for nutrients in nutrients_options[end-3:end], # empty slot, nitrate/ammonia, and N
        detritus in detritus_options,
        inorganic_carbon in inorganic_carbon_options,
        oxygen in oxygen_options

        biogeochemistry = NutrientsPlanktonDetritus(grid;
                                                    plankton = Abiotic(),
                                                    nutrients,
                                                    detritus,
                                                    inorganic_carbon,
                                                    oxygen)

        @info summary(biogeochemistry.underlying_biogeochemistry)

        model = NonhydrostaticModel(grid; 
                                    advection = nothing,
                                    biogeochemistry)

        time_step!(model, 1.0) 
        @test (CUDA.@allowscalar all([all(values .== 0) for values in values(model.tracers)]))

        check_conservations(model)
    end

    for plankton in (PhytoZoo, ImplicitProductivity),
        nutrients in nutrients_options[end-3:end], # empty slot, nitrate/ammonia, and N
        detritus in detritus_options,
        inorganic_carbon in inorganic_carbon_options,
        oxygen in oxygen_options

        nutrient_half_saturations = 
            NamedTuple(name=>half_sat for 
                        (name, half_sat) in pairs((phosphate = 1, iron = 0.001)) 
                        if !isnothing(getproperty(nutrients, name)))

        if nutrients.nitrogen isa NitrateAmmonia
            nutrient_half_saturations = merge(nutrient_half_saturations, (nitrate = 0.5, ammonia = 0.1))
        elseif nutrients.nitrogen isa OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels.SingleTracerNutrient
            nutrient_half_saturations = merge(nutrient_half_saturations, (; nitrate = 0.5))
        end

        biogeochemistry = NutrientsPlanktonDetritus(grid;
                                                    plankton = plankton(; nutrient_half_saturations),
                                                    nutrients,
                                                    detritus,
                                                    inorganic_carbon,
                                                    oxygen,
                                                    light_attenuation)

        @info summary(biogeochemistry.underlying_biogeochemistry)

        model = NonhydrostaticModel(grid; 
                                    advection = nothing,
                                    biogeochemistry)

        time_step!(model, 1.0) 
        @test (CUDA.@allowscalar all([all(values .== 0) for values in values(model.tracers)]))

        check_conservations(model)

        if !((detritus isa InstantRemineralisationDetritus)&(plankton == ImplicitProductivity))
            check_stepped(model)
        end
    end
end

@testset "Alkalinity sign" begin
    biogeochemistry = LOBSTER(grid; inorganic_carbon = CarbonateSystem(), plankton = PhytoZoo(rain_ratio = 0))

    model = NonhydrostaticModel(grid; advection = nothing, biogeochemistry)

    set!(model, NO₃ = 10, NH₄ = 0, P = 1, Z = 0.1, DOM = 0.1, sPOM = 0.1, bPOM = 0.1, DIC = 2000, Alk = 1900)

    NO₃₀, NH₄₀, Alk₀ = CUDA.@allowscalar model.tracers.NO₃[1, 1, 1], model.tracers.NH₄[1, 1, 1], model.tracers.Alk[1, 1, 1]
    time_step!(model, 1) 
    ΔNO₃, ΔNH₄, ΔAlk = CUDA.@allowscalar model.tracers.NO₃[1, 1, 1] - NO₃₀, model.tracers.NH₄[1, 1, 1] - NH₄₀, model.tracers.Alk[1, 1, 1] - Alk₀

    @test ΔNO₃ < 0 && ΔAlk > 0 # net nitrate uptake raises alkalinity
    @test isapprox(ΔAlk, -(1 + 1/16) * ΔNO₃ + (1 - 1/16) * ΔNH₄)
end

@testset "Explicit particle sinking" begin # not sure how essential this is but seems worth doing
    grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 2))

    biogeochemistry = NPZD(grid; 
                           detritus = Detritus(grid;
                                               sinking_speed=1,
                                               remineralisation_rate=0.0,
                                               open_bottom = true))

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, D = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.D[1, 1, 2], exp(-0.01*100), atol = 1e-3)
    end 

    biogeochemistry = NPZD(grid; 
                           detritus = Detritus(grid;
                                               sinking_speed=1,
                                               remineralisation_rate=0.0,
                                               open_bottom = false))

    # we scale to zero at the bottom which in a normal situation is fine but in this situation means the top cell is also not as set
    first_cell_flux_velocity = CUDA.@allowscalar biogeochemistry.underlying_biogeochemistry.detritus.sinking_speeds.w[1, 1, 2]

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, D = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.D[1, 1, 2], exp(0.01*100 * first_cell_flux_velocity), atol = 1e-3)
        @test Field(Integral(model.tracers.D))[1, 1, 1] ≈ 1
    end  

    biogeochemistry = NPZD(grid; 
                           detritus = DissolvedParticulate(grid;
                                                           sinking_speeds=(0.1, 1),
                                                           particulate_remineralisation_rate=(0.0, 0.0),
                                                           open_bottom = true))

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, sPOM = (x, y, z) -> z>-1, bPOM = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.sPOM[1, 1, 2], exp(-0.01*100*0.1), atol = 1e-3)
        @test isapprox(model.tracers.bPOM[1, 1, 2], exp(-0.01*100), atol = 1e-3)
    end 
end

using OceanBioME.Models.NutrientsPlanktonDetritusModels: SingleTracerNutrient

@testset "NutrientsPlanktonDetritus constructors and float types" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

    lobster = LOBSTER(grid)
    npzd = NPZD(grid)
    implicit = ImplicitBiology(grid)

    @test lobster isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test npzd isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test implicit isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}

    @test lobster.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:NitrateAmmonia, Nothing, Nothing, Nothing}, <:PhytoZoo, <:DissolvedParticulate{1, 2}}
    @test npzd.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:SingleTracerNutrient, Nothing, Nothing, Nothing}, <:PhytoZoo, <:Detritus}
    @test implicit.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:SingleTracerNutrient, <:SingleTracerNutrient, <:SingleTracerNutrient, Nothing}, <:ImplicitProductivity, <:DissolvedParticulate{1, 1}}

    grid = RectilinearGrid(architecture, Float32; size=(1, 1, 1), extent=(1, 1, 2))

    lobster = LOBSTER(grid; oxygen = Oxygen(Float32))
    npzd = NPZD(grid)
    implicit = ImplicitBiology(grid)

    @test lobster isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test npzd isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test implicit isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}

    @test lobster.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:NitrateAmmonia, Nothing, Nothing, Nothing}, <:PhytoZoo, <:DissolvedParticulate{1, 2}}
    @test npzd.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:SingleTracerNutrient, Nothing, Nothing, Nothing}, <:PhytoZoo, <:Detritus}
    @test implicit.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:SingleTracerNutrient, <:SingleTracerNutrient, <:SingleTracerNutrient, Nothing}, <:ImplicitProductivity, <:DissolvedParticulate{1, 1}}

    # light is converted
    par = lobster.light_attenuation
    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32
    @test par.phytoplankton_chlorophyll_ratio isa Float32

    # nutrients is converted
    @test lobster.underlying_biogeochemistry.nutrients.nitrogen |> ((::NitrateAmmonia{FT}) where FT) -> FT == Float32

    # plankton is converted
    @test lobster.underlying_biogeochemistry.plankton |> ((::PhytoZoo{<:Any, <:Any, <:Any, FT}) where FT) -> FT == Float32
    @test implicit.underlying_biogeochemistry.plankton |> ((::ImplicitProductivity{FT}) where FT) -> FT == Float32

    # detritus
    @test npzd.underlying_biogeochemistry.detritus |> ((::Detritus{FT}) where FT) -> FT == Float32
    @test lobster.underlying_biogeochemistry.detritus |> ((::DissolvedParticulate{<:Any, <:Any, <:Any, <:Any, FT}) where FT) -> FT == Float32

    # oxygen (carbonate system doesn't carry any numbers)
    @test lobster.underlying_biogeochemistry.oxygen |> ((::Oxygen{FT}) where FT) -> FT == Float32
end
