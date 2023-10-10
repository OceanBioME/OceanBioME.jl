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
using OceanBioME, Oceananigans.Units

const year = years = 365day
nothing #hide

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # specify the nominal depth of the box for the PAR profile
PAR(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay
nothing #hide

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid()), forcing = (; PAR))
model.Δt = 5minutes
model.stop_time = 5years

set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(model, save_interval = 100, save = SaveBoxModel("box.jld2"))

# ## Load the output

# Check the dependencies for loading and plotting the data are installed
# ```julia
# using Pkg
# pkg"add JLD2, CairoMakie"
# ```

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
using CairoMakie

fig = Figure(resolution = (800, 1600), fontsize = 24)

axs = []
for (idx, tracer) in enumerate(vars)
    push!(axs, Axis(fig[idx, 1], ylabel = "$tracer", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, timeseries[tracer], linewidth = 3)
end

fig
