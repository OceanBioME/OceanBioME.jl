# # Box model
# In this example we will setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration. This demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Plots, DiffEqBase, OrdinaryDiffEq"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME,  Plots, DiffEqBase, OrdinaryDiffEq 

# ## Define some convenient variables
day=days=60*60*24  # define the length of a day in seconds
year=years=day*365  # define the length of a year in days

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2

z=-10 # specify the nominal depth of the box for the PAR profile
PAR(t) = PAR⁰(t)*exp(z*0.2) # Modify the PAR based on the nominal depth and exponential decay 
params = merge(LOBSTER.defaults, (;PAR)) # Add the PAR to the list of parameters. The other parameters will be set to the defaults in the LOBSTER model

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = Setup.BoxModel(:LOBSTER, params, (NO₃=10.0, NH₄=0.1, P=0.1, Z=0.01, D=0.0, DD=0.0, Dᶜ=0.0, DDᶜ=0.0, DOM=0.0), 0.0, 2.0*year)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
solution = BoxModel.run(model) # call BoxModel to timestep the biogeochemical model

@info "Plotting the results..."
# ## Plot the results
values = vcat(transpose.(solution.u)...)

plts=[] # Create an empty array where we will store the sub-plots
# ## Loop over all tracers in the list, making a plot of the timeseries of each
for (i, tracer) in enumerate([:NO₃, :NH₄, :P, :Z, :DOM, :D, :DD])
    push!(plts, plot(solution.t/day, values[:, i], ylabel=tracer, xlabel="Day", legend=false))
end
# Add a plot of the PAR
push!(plts, plot(solution.t/day, PAR.(solution.t), xlabel="Day", ylabel="PAR",  legend=false))

plot(plts...) # display the combined plot

savefig("box.png") # save the plot to a file
