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
using OceanBioME,  Plots, DiffEqBase, OrdinaryDiffEq 
day=days=60*60*24  # define the length of a day in seconds
year=years=day*365  # define the length of a year in days

PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2

z=-10 # specify the depth for the light level
PAR(t) = PAR⁰(t)*exp(z*0.2) # Set the PAR
params = merge(LOBSTER.defaults, (;PAR))

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = Setup.BoxModel(:LOBSTER, params, (NO₃=10.0, NH₄=0.1, P=0.1, Z=0.01, D=0.0, DD=0.0, Dᶜ=0.0, DDᶜ=0.0, DOM=0.0), 0.0, 2.0*year)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
solution = BoxModel.run(model) # call BoxModel to timestep the biogeochemical model

@info "Plotting the results..."
# ## Plot the results
values = vcat(transpose.(solution.u)...)

plts=[]
for (i, tracer) in enumerate([:NO₃, :NH₄, :P, :Z, :DOM, :D, :DD])
    push!(plts, plot(solution.t/day, values[:, i], ylabel=tracer, xlabel="Day", legend=false))
end
push!(plts, plot(solution.t/day, PAR.(solution.t), xlabel="Day", ylabel="PAR",  legend=false))
plot(plts...)
savefig("examples/box.png")
