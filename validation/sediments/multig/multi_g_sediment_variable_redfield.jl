using Oceananigans, OceanBioME, Oceananigans.Units

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (16, ), extent = (10, ))

sediment_model =
    SimpleMultiGSediment(grid;
                         sinking_nitrogen = (:sPON, :bPON),
                         sinking_carbon = (:sPOC, :bPOC))

carbon_conservation = ScaleNegativeTracers((tracers = (:bPOC, :sPOC, :DOC, :DIC, :P, :Z),
                                            scalefactors = (1, 1, 1, 1, 6.56, 6.56)),
                                            grid)

biogeochemistry = LOBSTER(; grid, sediment_model, oxygen = true, variable_redfield = true, carbonates = true, scale_negatives = true, modifiers = carbon_conservation)

model = NonhydrostaticModel(grid; biogeochemistry, advection = WENO(order = 3))#, bounds = (0, 1)))

set!(model, sPOC = 6.56, bPOC = 6.56, sPON = 1, bPON = 1, O₂ = 500, NH₄ = 1, NO₃ = 10, DIC = 1500)

set!(sediment_model.fields.Ns, 10)
set!(sediment_model.fields.Nf, 1)

set!(sediment_model.fields.Cs, 10 * 6.56)
set!(sediment_model.fields.Cf, 1 * 6.56)

tracer_nitrogen = Field(Integral(model.tracers.sPON)) +
                  Field(Integral(model.tracers.bPON)) +
                  Field(Integral(model.tracers.P)) +
                  Field(Integral(model.tracers.Z)) +
                  Field(Integral(model.tracers.NO₃)) +
                  Field(Integral(model.tracers.NH₄)) +
                  Field(Integral(model.tracers.DON))

total_nitrogen = Field(tracer_nitrogen + sediment_model.fields.Nf + sediment_model.fields.Ns + sediment_model.fields.Nr)

tracer_carbon = Field(Integral(model.tracers.sPOC)) +
                Field(Integral(model.tracers.bPOC)) +
                Field(Integral(model.tracers.P)) * 6.56 +
                Field(Integral(model.tracers.Z)) * 6.56 +
                Field(Integral(model.tracers.DIC)) +
                Field(Integral(model.tracers.DOC))

total_carbon = Field(tracer_carbon + sediment_model.fields.Cf + sediment_model.fields.Cs + sediment_model.fields.Cr)

simulation = Simulation(model, Δt = 2minutes, stop_time = 4days)

compute!(total_nitrogen)
compute!(total_carbon)

@info total_nitrogen[1, 1, 1]
@info total_carbon[1, 1, 1]

simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers,
                                                 overwrite_existing = true,
                                                 filename = "sediment_tracers_redfield.jld2",
                                                 schedule = TimeInterval(30minutes))

simulation.output_writers[:sediment] = JLD2Writer(model, sediment_model.fields,
                                                  overwrite_existing = true,
                                                  filename = "sediment_redfield.jld2",
                                                  schedule = TimeInterval(30minutes))

prog(sim) = @info "$(prettytime(sim)) in $(prettytime(sim.run_wall_time))"

add_callback!(simulation, prog, IterationInterval(100))
#=
run!(simulation)

compute!(total_nitrogen)
compute!(total_carbon)

@info total_nitrogen[1, 1, 1]
@info total_carbon[1, 1, 1]
=#
