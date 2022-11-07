# # Simple active particle example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and active particles modelling the growth of sugar kelp. This demonstraits:
# - How to setup OceanBioME's biogeochemical models
# - How to setup light attenuation
# - How to add biologically active particles which interact with the biodeochemical model
# - How to include optional tracer sets (carbonate chemistry and oxygen)
# - How to visulise results

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, Printf, GLMakie2"
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

t_function(x, y, z, t) = 2.4*cos(t*2π/year + 50day) + 10
s_function(x, y, z, t) = 35.0

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# ## Kelp Particle setup
@info "Setting up kelp particles"
n_kelp = 5 # number of kelp fronds
z₀ = [-21:5:-1;]*1.0 # depth of kelp fronds

kelp_particles = SLatissima.setup(n_kelp, Lx/2, Ly/2, z₀, 
                                                    30.0, 0.01, 0.1, 57.5;
                                                    scalefactor = 500.0, 
                                                    T = t_function, S = s_function, urel = 0.2, 
                                                    optional_tracers = (:NH₄, :DIC, :DD, :DDᶜ, :OXY, :DOM))

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))
# ## Biogeochemical and Oceananigans model
# Here we instantiate the simplest possible version of the LOBSTER model which will return all of the information we then need to pass onto Oceananigans and set the initial conditions.
# > We are being a bit picky about the Oceananigans model setup (e.g. specifying the advection scheme) as this gives the best simple results but you may find differently.
bgc = Setup.Oceananigans(:LOBSTER, grid, params, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc)) 

@info "Setup BGC model"

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
set!(model, P=0.03, Z=0.03, D=0.0, DD=0.0, Dᶜ=0.0, DDᶜ=0.0, NO₃=11, NH₄=0.05, DOM=0.0, DIC=2200.0, ALK=2400.0, OXY=240.0)


# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=10minutes, stop_time=100days) 

simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "kelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day))
simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), filename = "$filename_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))

# ## Run!
# Finally we run the simulation
run!(simulation)
