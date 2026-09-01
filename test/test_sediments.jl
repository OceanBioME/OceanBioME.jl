include("dependencies_for_runtests.jl")

using OceanBioME.Models.SedimentModels: InstantRemineralisation, SimpleMultiG
using OceanBioME.Sediments: BiogeochemicalSediment

display_name(::NutrientsPlanktonDetritus) = "NutrientsPlanktonDetritus"
display_name(::BiogeochemicalSediment{<:SimpleMultiG}) = "Multi-G"
display_name(::BiogeochemicalSediment{<:OceanBioME.Models.SedimentModels.InstantRemineralisation}) = "Instant remineralisation"
display_name(::RectilinearGrid) = "Rectilinear grid"
display_name(::LatitudeLongitudeGrid) = "Latitude longitude grid"
display_name(::ImmersedBoundaryGrid) = "Immersed boundary grid"

function display_name(architecture, grid, sediment_model, biogeochemistry, timestepper)
    arch_name = typeof(architecture)
    sediment_name = display_name(sediment_model)
    bgc_name = display_name(biogeochemistry.underlying_biogeochemistry)
    grid_name = display_name(grid)

    @info "Testing sediment on $arch_name with $timestepper and $sediment_name on $bgc_name with $grid_name"

    return "$architecture, $timestepper, $sediment_name, $bgc_name, $grid_name"
end

set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Detritus}, model) = set!(model, D = 1)
set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:DissolvedParticulate{1, 2}}, model) = set!(model, sPOM = 1, bPOM = 1)
#=set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:DissolvedPar}, model) =
    set!(model, sPON = 1, bPON = 1, sPOC = 6.56, bPOC = 6.56)=#

sum_of_volume_integrals(biogeochemistry, tracers) = sum(map(f -> Field(Integral(f)), values(tracers)))
#=sum_of_volume_integrals(::NutrientsPlanktonDetritus{<:Any, <:Any, <:VariableRedfieldDetritus}, tracers) =
    sum([Field(Integral(f)) for (n, f) in pairs(tracers) if n in (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)])=#

sum_of_area_integrals(sediment, fields) = sum(map(f -> Field(Integral(f, dims = (1, 2))), values(fields)))
sum_of_area_integrals(::SimpleMultiG{Nothing}, fields) =
    sum([Field(Integral(f, dims = (1, 2))) for (n, f) in pairs(fields) if n in (:Nf, :Ns, :Nr)])

function test_sediment(grid, biogeochemistry, model_name, advection = WENO(order = 3, bounds = (0, 1)))
    method = quote
        return $(model_name)($(grid);
                             biogeochemistry = $(biogeochemistry),
                             buoyancy = nothing,
                             tracers = (),
                             $(ifelse(model_name == NonhydrostaticModel, :advection, :tracer_advection)) = $advection)
    end

    model = eval(method)

    sediment_model = biogeochemistry.sediment

    set_sinkers!(biogeochemistry.underlying_biogeochemistry, model)

    if isa(biogeochemistry.sediment.biogeochemistry, SimpleMultiG)
        set!(model, NO₃ = 10, NH₄ = 1, O₂ = 1000)
    end

    tracer_nitrogen = sum_of_volume_integrals(biogeochemistry.underlying_biogeochemistry, model.tracers)

    sediment_nitrogen = sum_of_area_integrals(biogeochemistry.sediment.biogeochemistry, sediment_model.fields)

    total_nitrogen = Field(tracer_nitrogen + sediment_nitrogen)

    compute!(total_nitrogen)

    initial_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

    for _ in 1:100
        time_step!(model, 1)
    end

    compute!(total_nitrogen)

    final_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

    # simple multi-G is only good to this precision, IR is fine to default
    @test isapprox(initial_total_nitrogen, final_total_nitrogen,rtol = 0.2e-6)

    @test all(interior(sediment_nitrogen) .!= 0)

    return model
end

bottom_height(x, y) = -1000 + 500 * exp(- (x^2 + y^2) / 250) # a perfect hill

rectilinear_grid = RectilinearGrid(
    architecture;
    size = (3, 3, 50),
    extent = (10, 10, 500)
)

latlon_grid = LatitudeLongitudeGrid(
    architecture;
    size = (3, 3, 16),
    latitude = (0, 10),
    longitude = (0, 10),
    z = (-500, 0)
)

immersed_latlon_grid = ImmersedBoundaryGrid(
    latlon_grid,
    GridFittedBottom(bottom_height)
)

grids = (rectilinear_grid, latlon_grid, immersed_latlon_grid)
sediment_timesteppers = (:QuasiAdamsBashforth2, :RungeKutta3)
models = (NonhydrostaticModel, HydrostaticFreeSurfaceModel) # I don't think we need to test on both models anymore
#=
@testset "Sediment integration" begin
    for grid in grids, timestepper in sediment_timesteppers
        npzd_ir = NutrientPhytoplanktonZooplanktonDetritus(;
            grid,
            sediment_model = InstantRemineralisationSediment(grid; timestepper)
        )

        lobster_ir = LOBSTER(;
            grid,
            sediment_model = InstantRemineralisationSediment(
                grid;
                sinking_tracers = (:sPOM, :bPOM),
                remineralisation_reciever = :NH₄,
                timestepper
            )
        )

        simple_lobster_multi_g = LOBSTER(;
            grid,
            sediment = SimpleMultiGSediment(grid),
            oxygen = true
        )

        full_lobster_multi_g = LOBSTER(;
            grid,
            sediment = SimpleMultiGSediment(
                grid;
                sinking_nitrogen = (:sPON, :bPON),
                sinking_carbon = (:sPOC, :bPOC)
            ),
            oxygen = true,
            carbonates = true,
            variable_redfield = true
        )

        bgcs = [npzd_ir, lobster_ir, simple_lobster_multi_g, full_lobster_multi_g]

        for model in models, biogeochemistry in bgcs
            nonhydrostatic = (model == NonhydrostaticModel)

            grid_is_immersed = grid isa ImmersedBoundaryGrid
            grid_is_latlon = grid isa LatitudeLongitudeGrid

            # Skip incompatible combinations
            if nonhydrostatic && (grid_is_immersed || grid_is_latlon)
                continue
            end

            test_name = display_name(architecture, grid, biogeochemistry.sediment, biogeochemistry, timestepper)

            @testset "$(test_name)" begin
                test_sediment(grid, biogeochemistry, model)
            end
        end
    end
end
=#

#####
##### The sedimentary denitrification cap needs the floor cell's thickness, and the floor is not always
##### k = 1. Build the same physical floor cell two ways — a rectilinear grid cut to the sea floor
##### (k_bottom = 1, the shape every MARBL gate has used) and a full-depth immersed grid whose bathymetry
##### puts the floor at the same depth (k_bottom > 1) — and require the tendencies to agree. The vertical
##### spacing is stretched *across the floor*, so a thickness taken from the face-centred spacing (which
##### straddles the floor cell and the one below it) differs from the cell's own and the test fails.
#####

using OceanBioME.Models.SedimentModels: BurialDenitrification, bottom_Δz

sediment_tracked_fields(grid; POC, O₂, NO₃) =
    NamedTuple{(:POC, :POP, :bSi, :CaCO₃, :PFe, :O₂, :NO₃, :Ω, :O₂_consumption_scale)}(
        map((POC, POC / 117, POC / 15, POC, POC / 3e4, O₂, NO₃, 1, 1)) do value
            field = Field{Center, Center, Nothing}(grid)
            set!(field, value)
            field
        end)

@testset "Denitrification clamp on an immersed floor" begin
    # stretched so that the floor cell and the cell below it have different thicknesses
    z_faces = [-500, -400, -200, -150, -100, -70, -45, -25, -10, 0]
    k_floor = 3    # the cell spanning -200 to -150
    z_floor = z_faces[k_floor]

    full_grid = RectilinearGrid(architecture; size = (1, 1, length(z_faces) - 1),
                                x = (0, 1), y = (0, 1), z = z_faces,
                                topology = (Bounded, Bounded, Bounded))

    cut_grid = RectilinearGrid(architecture; size = (1, 1, length(z_faces) - k_floor),
                               x = (0, 1), y = (0, 1), z = z_faces[k_floor:end],
                               topology = (Bounded, Bounded, Bounded))

    # a mid-cell floor, where the cell is cut to a thickness no face spacing knows about
    z_partial = (z_faces[k_floor] + z_faces[k_floor + 1]) / 2

    partial_cut_grid = RectilinearGrid(architecture; size = (1, 1, length(z_faces) - k_floor),
                                       x = (0, 1), y = (0, 1),
                                       z = [z_partial, z_faces[(k_floor + 1):end]...],
                                       topology = (Bounded, Bounded, Bounded))

    immersed_grids = (GridFittedBottom(z_floor) => cut_grid,
                      PartialCellBottom(z_partial) => partial_cut_grid)

    sediment = BurialDenitrification()

    # suboxic and nitrate poor, so the cap actually binds
    POC, O₂, NO₃ = 5e-4, 1.0, 0.4

    tendency(grid) =
        sediment(1, 1, grid, Val(:denitrified_N), nothing, nothing,
                 sediment_tracked_fields(grid; POC, O₂, NO₃))

    # a cap of exactly one would make every comparison below vacuous
    clamp = OceanBioME.Models.SedimentModels.denitrification_clamp(sediment, POC, O₂, NO₃,
                                                                  bottom_Δz(cut_grid, 1, 1))

    @test 0 < clamp < 1

    for (bottom, equivalent_grid) in immersed_grids
        immersed_grid = ImmersedBoundaryGrid(full_grid, bottom)

        expected = CUDA.@allowscalar bottom_Δz(equivalent_grid, 1, 1)

        @test (CUDA.@allowscalar bottom_Δz(immersed_grid, 1, 1)) ≈ expected

        # the thickness must be the one the coupling divides the returned area flux by
        @test expected ≈ CUDA.@allowscalar Oceananigans.Operators.Δzᶜᶜᶜ(1, 1, 1, equivalent_grid)

        @test (CUDA.@allowscalar tendency(immersed_grid)) ≈ CUDA.@allowscalar tendency(equivalent_grid)
    end
end

#####
##### The implicit ballast sweep starts each column at its own sea floor. That index is built by
##### `floor_index_field`, which had never run on an immersed grid, so check both that it does and that the
##### sweep it feeds completes a step — the index is a loop bound (`for k in Nz:-1:k_bottom`) and is
##### compared to `k`, so it has to come back as an integer.
#####

using OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels.MultiElement: floor_index_field

@testset "Implicit ballast floor index on an immersed grid" begin
    z_faces = [-500, -400, -200, -150, -100, -70, -45, -25, -10, 0]

    grid = RectilinearGrid(architecture; size = (4, 4, length(z_faces) - 1), halo = (3, 3, 3),
                           x = (0, 1), y = (0, 1), z = z_faces,
                           topology = (Bounded, Bounded, Bounded))

    # a two step sea floor, so the columns start the sweep at different levels and a k = 1 assumption shows
    stepped_bottom(x, y) = ifelse(x < 0.5, -200, -100)

    bottoms = (GridFittedBottom(stepped_bottom),
               PartialCellBottom((x, y) -> stepped_bottom(x, y) + 25))

    for bottom in bottoms
        immersed_grid = ImmersedBoundaryGrid(grid, bottom)

        floor_indices = floor_index_field(immersed_grid)

        @test eltype(floor_indices) <: Integer

        # -200 is the base of the third cell, -100 the base of the fifth
        @test Array(interior(floor_indices, :, 1, 1)) == [3, 3, 5, 5]
    end

    # ...and the sweep those indices drive runs to completion on an immersed column
    immersed_grid = ImmersedBoundaryGrid(grid, GridFittedBottom(stepped_bottom))

    # PAR is prescribed, so the sweep is what is under test rather than the light model
    biogeochemistry = MARBL(immersed_grid;
                            light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(50.0)))

    model = HydrostaticFreeSurfaceModel(immersed_grid; biogeochemistry,
                                        tracer_advection = WENO(order = 3),
                                        buoyancy = nothing, tracers = (:T, :S))

    set!(model, T = 10, S = 35, DIC = 2000, Alk = 2300, NO₃ = 10, NH₄ = 0.1, PO₄ = 1,
         Si = 10, Fe = 1e-3, Lig = 1e-3, O₂ = 200,
         spC = 1, spChl = 0.03, spP = 1e-2, spFe = 5e-5, spCaCO₃ = 1e-2,
         diatC = 1, diatChl = 0.03, diatP = 1e-2, diatFe = 5e-5, diatSi = 0.1,
         diazC = 0.1, diazChl = 3e-3, diazP = 1e-3, diazFe = 5e-6,
         zooC = 1, DOC = 10, DON = 1, DOP = 0.1, DOCr = 10, DONr = 1, DOPr = 0.1)

    time_step!(model, 100)

    @test model.clock.iteration == 1

    detritus = biogeochemistry.underlying_biogeochemistry.detritus

    floor_flux = Array(interior(detritus.floor_flux.POC, :, 1, 1))
    remineralisation = Array(interior(detritus.remineralisation.POC, 1, 1, :))

    @test all(isfinite, floor_flux)
    @test all(isfinite, remineralisation)

    # particulate organic carbon is produced, so it reaches the floor of every column
    @test all(floor_flux .> 0)

    # the sweep starts at the floor, so it never writes to the immersed cells below it
    @test all(remineralisation[1:2] .== 0)
    @test any(remineralisation[3:end] .!= 0)

    @info "floor indices $(Array(interior(detritus.floor_indices, :, 1, 1))) " *
          "($(eltype(detritus.floor_indices)))\n  floor POC flux $(floor_flux)" *
          "\n  POC remineralisation (column 1) $(remineralisation)"
end
