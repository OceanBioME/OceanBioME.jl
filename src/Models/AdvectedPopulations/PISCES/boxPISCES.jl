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
include("PISCES.jl")

const year = years = 365day
nothing #hide

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # specify the nominal depth of the box for the PAR profile
PAR(t) = 3000# Modify the PAR based on the nominal depth and exponential decay
T(t) = 14.0

#PAR¹(t) = 100.0
#PAR²(t) = 100.0
#PAR³(t) = 100.0
#zₘₓₗ(t) = 100.0
#zₑᵤ(t) = 50.0

nothing #hide

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(biogeochemistry = PISCES(grid = BoxModelGrid, light_attenuation_model = nothing), forcing = (; PAR, T)) #, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ

set!(model,NO₃ = 4.0, NH₄ = 0.1, P = .426, D = .426, Z = 0.426, M = 0.426, Pᶠᵉ = 3.5, Dᶠᵉ = 3.5, Pᶜʰˡ = 1.60, Dᶜʰˡ = 1.60, Dˢⁱ = 0.525, Fe = 0.8241, O₂ = 264.0, Si = 4.557, Alk = 2360.0, PO₄ = 0.8114, DIC = 2000.0, CaCO₃ = 1.0, PAR¹ = 1000.0, PAR² = 1000.0,  PAR³ = 1000.0, zₘₓₗ = 100.0, zₑᵤ = 50.0)

simulation = Simulation(model; Δt = 0.5, stop_time = 12hours)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(300), overwrite_existing = true)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

simulation.callbacks[:progress] = Callback(prog, IterationInterval(10))

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(simulation)

# ## Load the output

times = FieldTimeSeries("box.jld2", "P").times

timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("box.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 7200), fontsize = 24)

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end

fig
