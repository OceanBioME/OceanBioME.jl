# # [Box model](@id box_example)
# In this example we setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME, Oceananigans, Oceananigans.Units

using Oceananigans.Fields: FunctionField, ConstantField

using OceanBioME.Models.PISCESModel: SimpleIron, NitrateAmmonia

const year = years = 365day

grid = RectilinearGrid( topology = (Flat, Flat, Flat), size = (), z = -10)
clock = Clock(time = zero(grid))

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
@inline surface_PAR(t) = 300 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 10
@inline temp(t) = 2.4 * (1 + cos(t * 2π / year + 50days)) + 8

@inline PAR_component(t) = surface_PAR(t) * exp(-10/2) / 3

PAR₁ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)
PAR₂ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)
PAR₃ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)

PAR = PAR₁ + PAR₂ + PAR₃

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

mixed_layer_depth = ConstantField(-100)
euphotic_depth = ConstantField(-100)
mean_mixed_layer_vertical_diffusivity = ConstantField(1e-2)
mean_mixed_layer_light = PAR
silicate_climatology = ConstantField(0) # turn off the contribution from enhanced requirments

biogeochemistry = PISCES(; grid, 
                           sinking_speeds = (POC = 0, GOC = 0), 
                           light_attenuation, 
                           mixed_layer_depth, 
                           euphotic_depth,
                           silicate_climatology,
                           mean_mixed_layer_light,
                           mean_mixed_layer_vertical_diffusivity,
                           iron = SimpleIron(excess_scavenging_enhancement = 0.0),
                           nitrogen = NitrateAmmonia(maximum_fixation_rate = 0.0))

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(; grid, biogeochemistry, clock, prescribed_tracers = (; T = temp))

set!(model, P = 0.1, PChl = 0.025, PFe = 0.005,
            D = 0.01, DChl = 0.003, DFe = 0.0006, DSi = 0.004,
            Z = 0.06, M = 0.5,
            DOC = 4, 
            POC = 5.4, SFe = 0.34, 
            GOC = 8.2, BFe = 0.5, PSi = 0.04, CaCO₃ = 10^-10,
            NO₃ = 10, NH₄ = 0.1, PO₄ = 5.0, Fe = 0.6, Si = 8.6,
            DIC = 2205, Alk = 2560, O₂ = 317, S = 35)

            
simulation = Simulation(model; Δt = 40minutes, stop_time = 4years)

simulation.output_writers[:fields] = JLD2Writer(model, model.fields; filename = "box.jld2", schedule = TimeInterval(0.5day), overwrite_existing = true)

PAR_field = Field(biogeochemistry.light_attenuation.fields[1])
simulation.output_writers[:par] = JLD2Writer(model, (; PAR = PAR_field); filename = "box_light.jld2", schedule = TimeInterval(0.5day), overwrite_existing = true)


prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

simulation.callbacks[:progress] = Callback(prog, TimeInterval(182days))

@info "Running the model..."
run!(simulation)

# load and plot results
timeseries = FieldDataset("box.jld2")

times = timeseries.fields["P"].times

PAR_timeseries = FieldTimeSeries("box_light.jld2", "PAR")

using CairoMakie

fig = Figure(size = (2400, 3600), fontsize = 24)

axs = []

n_start = 1

for name in Oceananigans.Biogeochemistry.required_biogeochemical_tracers(biogeochemistry)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/4), Int(idx%4)], ylabel = "$name", xlabel = "years", xticks=(0:40)))
    lines!(axs[end], times[n_start:end] / year, timeseries["$name"][n_start:end], linewidth = 3)
end

push!(axs, Axis(fig[6, 2], ylabel = "PAR", xlabel = "years", xticks=(0:40)))
lines!(axs[end], times[n_start:end]/year, PAR_timeseries[n_start:end], linewidth = 3)

fig
