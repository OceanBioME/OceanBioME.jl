# # One dimensional column example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstraits:
# - How to setup OceanBioME's biogeochemical models
# - How to setup light attenuation
# - How to visulise results

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, Printf, Plots, GLMakie, NetCDF, JLD2"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans,Printf
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
params = LOBSTER.defaults  

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((t-200days)/50days)^2))) .+ 2

H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
grid = RectilinearGrid(size=(1, 1, 50), extent=(20, 20, 200)) 
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# ## Biogeochemical and Oceananigans model
# Here we instantiate the simplest possible version of the LOBSTER model which will return all of the information we then need to pass onto Oceananigans and set the initial conditions.
# > We are being a bit picky about the Oceananigans model setup (e.g. specifying the advection scheme) as this gives the best simple results but you may find differently.
bgc = Setup.Oceananigans(:LOBSTER, grid, params) 

model = NonhydrostaticModel(
                                                advection = WENO(;grid),
                                                timestepper = :RungeKutta3,
                                                grid = grid,
                                                tracers = bgc.tracers,
                                                closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                                                forcing = bgc.forcing,
                                                boundary_conditions = bgc.boundary_conditions,
                                                auxiliary_fields = (; PAR)
)
set!(model, P=0.03, Z=0.03, D=0.0, DD=0.0, NO₃=11, NH₄=0.05, DOM=0.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=10minutes, stop_time=100days) 

simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)));

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day))
simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))

# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visulise the results

P = FieldTimeSeries("$filename.jld2", "P")
NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
Z = FieldTimeSeries("$filename.jld2", "Z")
D = FieldTimeSeries("$filename.jld2", "D") .+ FieldTimeSeries("$filename.jld2", "DD")
x, y, z = nodes(P)
times = P.times

using GLMakie
f=Figure(backgroundcolor=RGBf(1, 1, 1), fontsize=30)

axP = Axis(f[1, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Phytoplankton concentration (mmol N/m³)")
hmP = GLMakie.heatmap!(times./days, z[35:50], P[1, 1, 35:50, 1:101]', interpolate=true, colormap=:batlow)
cbP = Colorbar(f[1, 3], hmP)

axNO₃ = Axis(f[1, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Nitrate concentration (mmol N/m³)")
hmNO₃ = GLMakie.heatmap!(times./days, z[35:50], NO₃[1, 1, 35:50, 1:101]', interpolate=true, colormap=:batlow)
cbNO₃ = Colorbar(f[1, 6], hmNO₃)

axZ = Axis(f[2, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Zooplankton concentration (mmol N/m³)")
hmZ = GLMakie.heatmap!(times./days, z[35:50], Z[1, 1, 35:50, 1:101]', interpolate=true, colormap=:batlow)
cbZ = Colorbar(f[2, 3], hmZ)

axD = Axis(f[2, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Detritus concentration (mmol N/m³)")
hmD = GLMakie.heatmap!(times./days, z[35:50], D[1, 1, 35:50, 1:101]', interpolate=true, colormap=:batlow)
cbD = Colorbar(f[2, 6], hmD)