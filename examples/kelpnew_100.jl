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
using Interpolations
params = LOBSTER.defaults  

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

# PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((t-200days)/50days)^2))) .+ 2
# H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
# fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
# MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
# κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 
########################################################## temperature
t_temperature_node=[0.,66,95,240,364]
temperature_idealize=[9.0,8.05,8.05,13.65,9.0]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) 
########################################################## 
########################################################## salinity
s_function(x, y, z, t) = 35.0
########################################################## 
########################################################## PAR
t_par_node=[0.,30,120,200,330,364]
par_idealize=[3,3,90,90,3,3]
PAR_itp = LinearInterpolation((t_par_node)days, par_idealize)
PAR⁰(t) = PAR_itp(mod(t, 364days)) 
########################################################## 
########################################################## diffusivity 
t_mldplus_node=[0.,55,85,100,300,364]
mldplus_idealize=[280,420,420,40,40,280]
mld_itp = LinearInterpolation((t_mldplus_node)days, mldplus_idealize)  #in seconds 
κₜ(x, y, z, t) = 8e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4;
##########################################################

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
duration= 6years

Lx, Ly, Lz = 20, 20, 600
Nx, Ny, Nz = 1, 1, 60
grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz)) 
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# ## Kelp Particle setup
@info "Setting up kelp particles"
n_kelp = 10 # number of kelp fronds  100
z₀ = [-95:10:-1;]*1.0 # depth of kelp fronds   [-100:-1;]*1.0 
kelp_particles = SLatissima.setup(n_kelp, Lx/2, Ly/2, z₀, 
                                                    0.1, 0.01, 0.1, 57.5;
                                                    scalefactor = 100.0, 
                                                    T = t_function, S = s_function, urel = 0.2, 
                                                    optional_tracers = (:NH₄, :DIC, :DD, :DDᶜ, :OXY, :DOM))

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function), parameters=merge(Boundaries.defaults.airseaflux, (conc_air = (CO₂ = 400,), )))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))
# ## Biogeochemical and Oceananigans model
# Here we instantiate the simplest possible version of the LOBSTER model which will return all of the information we then need to pass onto Oceananigans and set the initial conditions.
# > We are being a bit picky about the Oceananigans model setup (e.g. specifying the advection scheme) as this gives the best simple results but you may find differently.
bgc = Setup.Oceananigans(:LOBSTER, grid, params, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), open_bottom=true) 

@info "Setup BGC model"

model = NonhydrostaticModel(
                                                advection = WENO(;grid),
                                                timestepper = :RungeKutta3,
                                                grid = grid,                  # 新的里面没有coriolis和buoyancy
                                                tracers = bgc.tracers,
                                                closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                                                forcing = bgc.forcing,
                                                boundary_conditions = bgc.boundary_conditions,
                                                auxiliary_fields = (; PAR),
                                                particles = kelp_particles
)
Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002           #in mmolN m^-3
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008           #in mmolN m^-3
Dᵢ(x, y, z)=0                                                #in mmolN m^-3
DDᵢ(x, y, z)=0                                               #in mmolN m^-3
Dᶜᵢ(x, y, z)=0                                                #in mmolC m^-3
DDᶜᵢ(x, y, z)=0                                               #in mmolC m^-3
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                #in mmolN m^-3
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in mmolN m^-3
DOMᵢ(x, y, z)= 0                                             #in mmolN m^-3
DICᵢ(x, y, z)= 2200                                          #in mmolC m^-3
ALKᵢ(x, y, z)= 2400                                          #in mmolN m^-3
OXYᵢ(x, y, z) = 240                                          #in mmolO m^-3
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, Dᶜ=Dᶜᵢ, DDᶜ=DDᶜᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, DIC=DICᵢ, ALK=ALKᵢ, OXY=OXYᵢ)  #没有速度了


# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=2.5minutes, stop_time=duration) 

# couple the particles and tracers
simulation.callbacks[:couple_particles] = Callback(OceanBioME.Particles.infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "kelpnew_100"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))
simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.1, c_adv=0.5, c_diff=0.5, w = 200/day, relaxation=0.75), TimeStepCallsite()) #没有kelp 可以不用c_forcing=0.5
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([i*years for i=1:12]), prefix=filename) #prefix="kelp_01"
# ## Run!
# Finally we run the simulation
run!(simulation)
