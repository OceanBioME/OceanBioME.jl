using Oceananigans, OceanBioME, Oceananigans.Units

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (16, ), extent = (10, ))

sediment_model = OceanBioME.Models.SedimentModels.SimpleMultiGSediment(grid;);

biogeochemistry = LOBSTER(; grid, sediment_model, oxygen = true, scale_negatives = true)

model = NonhydrostaticModel(; grid, biogeochemistry, advection = WENO(order = 3, bounds = (0, 1)))

set!(model, sPOM = 1, bPOM = 1, O₂ = 500, NH₄ = 1, NO₃ = 10)
set!(sediment_model.fields.Ns, 10)
set!(sediment_model.fields.Nf, 1)

tracer_nitrogen = Field(Integral(model.tracers.sPOM)) + 
                  Field(Integral(model.tracers.bPOM)) + 
                  Field(Integral(model.tracers.P)) + 
                  Field(Integral(model.tracers.Z)) + 
                  Field(Integral(model.tracers.NO₃)) + 
                  Field(Integral(model.tracers.NH₄)) + 
                  Field(Integral(model.tracers.DOM))

total_nitrogen = Field(tracer_nitrogen + sediment_model.fields.Nf + sediment_model.fields.Ns + sediment_model.fields.Nr)

compute!(total_nitrogen)

simulation = Simulation(model, Δt = 2minutes, stop_time = 4days)

@info total_nitrogen[1, 1, 1]

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       overwrite_existing = true, 
                                                       filename = "sediment_tracers.jld2", 
                                                       schedule = TimeInterval(30minutes))

simulation.output_writers[:sediment] = JLD2OutputWriter(model, sediment_model.fields,
                                                        overwrite_existing = true, 
                                                        filename = "sediment.jld2", 
                                                        schedule = TimeInterval(30minutes))

prog(sim) = @info "$(prettytime(sim)) in $(prettytime(sim.run_wall_time))"

add_callback!(simulation, prog, IterationInterval(100))
#=
run!(simulation)

compute!(total_nitrogen)

@info total_nitrogen[1, 1, 1]
=#