# # Box model
# In this example we will setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration. This demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, CairoMakie, DiffEqBase, OrdinaryDiffEq"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME

minute=minutes=60
hour=hours=60*minutes
day=days=hours*24  # define the length of a day in seconds
year=years=day*365  # define the length of a year in days

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2

z=-10# specify the nominal depth of the box for the PAR profile
PAR(t) = PAR⁰(t)*exp(z*0.2) # Modify the PAR based on the nominal depth and exponential decay 

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid()), forcing = (; PAR))
model.Δt = 5minutes
model.stop_time = 2years

set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(model, save_interval = 100, save = SaveBoxModel("box.jld2"))

@info "Plotting the results..."
# ## Plot the results
using JLD2, CairoMakie
vars = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :PAR)
file = jldopen("box.jld2")
times = keys(file["values"])
timeseries = NamedTuple{vars}(ntuple(t -> zeros(length(times)), length(vars)))

for (idx, time) in enumerate(times)
    values = file["values/$time"]
    for tracer in vars
        getproperty(timeseries, tracer)[idx] = values[tracer]
    end
end

close(file)

fig = Figure(resolution = (1600, 1000))

plt_times = parse.(Float64, times)./day

axs = []

for (idx, tracer) in enumerate(vars)
    push!(axs, Axis(fig[floor(Int, (idx - 1)/4) + 1, (idx - 1) % 4 + 1], ylabel="$tracer", xlabel="Day"))
    lines!(axs[end], plt_times, timeseries[tracer])
end
save("box.png", fig)

# ![Results](box.png)
