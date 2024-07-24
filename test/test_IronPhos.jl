include("dependencies_for_runtests.jl")

using OceanBioME.Sediments: SimpleMultiG, InstantRemineralisation, IronPhosphate
using Oceananigans.Units
using Oceananigans
using OceanBioME.Sediments: sediment_tracers, sediment_fields
using Oceananigans: Field
using Oceananigans.Fields: TracerFields

using Oceananigans.Operators: volume, Azᶠᶜᶜ

using OceanBioME.LOBSTERModel: VariableRedfieldLobster

function intercept_tracer_tendencies!(model, intercepted_tendencies)
    for (name, field) in enumerate(intercepted_tendencies)
        field .= Array(interior(model.timestepper.Gⁿ[name + 3]))
    end
end

function set_defaults!(sediment::IronPhosphate)
    set!(sediment.fields.O₂, 1e-6)
    set!(sediment.fields.NH₄, 100e-6)
    set!(sediment.fields.NO₃, 6e-6)
    set!(sediment.fields.NO₂, 0.5e-6)
    set!(sediment.fields.N₂, 0)
    set!(sediment.fields.TPO₄, 1.5e-6)
    set!(sediment.fields.FeOHP, 410e-6)
    set!(sediment.fields.Feᴵᴵ, 0.1e-6)
    set!(sediment.fields.FeS₂, 0)
    set!(sediment.fields.SO₄, 0.5e-6)
    set!(sediment.fields.TH₂S, 0.1e-6)
    set!(sediment.fields.CH₄, 10e-9)
    set!(sediment.fields.TCO₂, 0)
    set!(sediment.fields.Gi, 0)
end 

set_defaults!(::VariableRedfieldLobster, model) =
    set!(model, P = 0.4686, Z = 0.5363, 
                NO₃ = 2.3103, NH₄ = 0.0010, 
                DIC = 2106.9, Alk = 2408.9, 
                O₂ = 258.92, 
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781)

#total_nitrogen(sed::SimpleMultiG) = sum(sed.fields.N_fast) + 
#                                    sum(sed.fields.N_slow) + 
#                                    sum(sed.fields.N_ref)

#total_nitrogen(::VariableRedfieldLobster, model) = sum(model.tracers.NO₃) +
#                                                     sum(model.tracers.NH₄) +
#                                                     sum(model.tracers.P) +
#                                                     sum(model.tracers.Z) +
#                                                     sum(model.tracers.DON) +
#                                                     sum(model.tracers.sPON) +
#                                                     sum(model.tracers.bPON)

bottom_height(x, y) = -1000 + 500 * exp(- (x^2 + y^2) / 250) # a perfect hill

grid = RectilinearGrid(architecture; size=(3, 3, 50), extent=(10, 10, 500))
sediment_model = IronPhosphate(; grid)
biogeochemistry = LOBSTER(; grid,
                                    carbonates = true, 
                                    oxygen = true, 
                                    variable_redfield = true, 
                                    sediment_model)

model = NonhydrostaticModel(; grid, 
            biogeochemistry, 
            closure = nothing,
            timestepper = :RungeKutta3,
            buoyancy = nothing)

set_defaults!(model.biogeochemistry.sediment)

set_defaults!(biogeochemistry.underlying_biogeochemistry, model)

sim_length = 5days

simulation = Simulation(model, Δt = 50, stop_time = sim_length)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.biogeochemistry.sediment.fields,
    filename = "temp_plotting_data.jld2",
    schedule = TimeInterval(24minute),
    overwrite_existing = true)

intercepted_tendencies = Tuple(Array(interior(field)) for field in values(TracerFields(keys(model.tracers), grid)))

simulation.callbacks[:intercept_tendencies] = Callback(intercept_tracer_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

run!(simulation)

var_name_example = keys(model.biogeochemistry.sediment.fields)[1]
times = FieldTimeSeries("temp_plotting_data.jld2", "$var_name_example").times

timeseries = NamedTuple{keys(model.biogeochemistry.sediment.fields)}(FieldTimeSeries("temp_plotting_data.jld2", "$field")[1, 1, 1, :] for field in keys(model.biogeochemistry.sediment.fields))

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 3600), fontsize = 24)

tick_location_seconds = range(0, times[end]; length=5)
tick_location_days = tick_location_seconds / (24 * 60 * 60)


axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Days", xticks=(collect(tick_location_seconds), string.(collect(tick_location_days)))))
    lines!(axs[end], times, tracer, linewidth = 3)
end

display(fig)



@info "Success!"