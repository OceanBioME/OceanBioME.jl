# # Box model
# In this example we will setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration. This demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, CairoMakie, JLD2"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME

const minute = minutes = 60
const hour = hours = 60 * minutes
const day = days = hours * 24  # define the length of a day in seconds
const year = years = day * 365 # define the length of a year in days

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # specify the nominal depth of the box for the PAR profile
PAR(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay 

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid()), forcing = (; PAR))
model.Δt = 5minutes
model.stop_time = 5years

set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(model, save_interval = 100, save = SaveBoxModel("box.jld2"))

# ## Load the output
@info "Loading output..."

using JLD2

vars = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :PAR)
file = jldopen("box.jld2")
times = parse.(Float64, keys(file["values"]))

timeseries = NamedTuple{vars}(ntuple(t -> zeros(length(times)), length(vars)))

for (idx, time) in enumerate(times)
    values = file["values/$time"]
    for tracer in vars
        getproperty(timeseries, tracer)[idx] = values[tracer]
    end
end

close(file)

# ## And plot
@info "Plotting the results..."

using CairoMakie

fig = Figure(resolution = (800, 1600), fontsize = 24)

axs = []
for (idx, tracer) in enumerate(vars)
    push!(axs, Axis(fig[idx, 1], ylabel = "$tracer", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, timeseries[tracer], linewidth = 3)
end

fig
