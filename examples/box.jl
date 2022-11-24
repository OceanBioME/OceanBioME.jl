# # Box model
# In this example we will setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration. This demonstraits:
# - How to setup OceanBioME's biogeochemical models as box models

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Plots, DiffEqBase, OrdinaryDiffEq"
# ```

# ## Model setup
# We load the packages and setup the initial and forcing conditions
using OceanBioME

import OceanBioME.BoxModels: update_boxmodel_state!

minute=minutes=60
hour=hours=60*minutes
day=days=hours*24  # define the length of a day in seconds
year=years=day*365  # define the length of a year in days

PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2

z=-10 # specify the depth for the light level
PAR(t) = PAR⁰(t)*exp(z*0.2) # Set the PAR

# Clumsily setup how we will modify the PAR in the model
function update_boxmodel_state!(model::BoxModel{<:LOBSTER, <:Any, <:Any, <:Any, <:Any, <:Any})
    getproperty(model.values, :PAR) .= PAR(model.clock.time)
end

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid()))
model.Δt = 5minutes
model.stop_time = 2years

set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(model, save_interval = 100, save = SaveBoxModel("box.jld2"))

@info "Plotting the results..."
# ## Plot the results
using JLD2, Plots
vars =(:NO₃, :NH₄, :P, :Z, :DOM, :D, :DD, :PAR)
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

plts=[]
for tracer in vars
    push!(plts, plot(parse.(Float64, times)./day, timeseries[tracer], ylabel=tracer, xlabel="Day", legend=false))
end
plot(plts...)
savefig("examples/box.png")
