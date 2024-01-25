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
# pkg"add OceanBioME, Oceananigans, Printf, GLMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans,Printf
using Oceananigans.Fields: FunctionField, ConstantField 
# using Oceananigans.Units: second, minute, minutes, hour, hours, day, days#, year, years
using Oceananigans.Units #####changed
using Interpolations
using JLD2
using Oceananigans.Architectures: arch_array  #####changed
const year = years = 365days                  #####changed

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

##PAR⁰(x, y, t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2
##H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
##fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
##MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
##κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 
##t_function(x, y, z, t) = 2.4*cos(t*2π/year + 50day) + 10
##s_function(x, y, z, t) = 35.0

########################################################## temperature
t_temperature_node=[0.,66,95,240,364]
temperature_idealize=[9.0,8.05,8.05,13.65,9.0]
temperature_itp = LinearInterpolation((t_temperature_node)days, temperature_idealize) 
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) 
########################################################## 
########################################################## salinity
s_function(x, y, z, t) = 35.16
########################################################## 
########################################################## PAR
t_par_node=[0.,30,120,200,330,364]
par_idealize=[3,3,90,90,3,3]
PAR_itp = LinearInterpolation((t_par_node)days, par_idealize)
PAR⁰(x, y, t) = PAR_itp(mod(t, 364days)) 
########################################################## 
########################################################## diffusivity 
t_mldplus_node=[0.,55,85,100,300,364]
mldplus_idealize=[280,420,420,40,40,280]
mld_itp = LinearInterpolation((t_mldplus_node)days, mldplus_idealize)  #in seconds 
κₜ(x, y, z, t) = 8e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4;
##########################################################
# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
##Lx, Ly = 20, 20
##grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 
##PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  
architecture = CPU()
duration= 2years
Lx, Ly, Lz = 20, 20, 600
Nx, Ny, Nz = 1, 1, 100
grid = RectilinearGrid(architecture, size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz))   ####change
#PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(t_function, grid; clock)
# S = ConstantField(35.16)  # has error: LoadError: type ConstantField has no field grid
S = FunctionField{Center, Center, Center}(s_function, grid; clock)

# ## Kelp Particle setup
@info "Setting up kelp particles"
n = 100 # number of kelp fronds
z₀ = [-100:-1;]*1.0 # depth of kelp fronds   z₀ = [-100:-1;]*1.0   [-21:5:-1;]*1.0

# kelp_particles = SLatissima.setup(;n, 
#                                   x₀ = Lx/2, y₀ = Ly/2, z₀, 
#                                   A₀ = 0.1, N₀ = 0.01, C₀ = 0.1, 
#                                   latitude = 57.5,
#                                   scalefactor = 450.0, 
#                                   T = t_function, S = s_function, urel = 0.2, 
#                                   optional_tracers = (:NH₄, :DIC, :bPON, :bPOC, :O₂, :DON, :DOC))
kelp_particles = SLatissima(; architecture, 
                                  x = arch_array(architecture, ones(n) * Lx / 2), 
                                  y = arch_array(architecture, ones(n) * Ly / 2), 
                                  z = arch_array(architecture, z₀), 
                                  A = arch_array(architecture, ones(n) * 0.1),
                                  N = arch_array(architecture, ones(n) * 0.01),
                                  C = arch_array(architecture, ones(n) * 0.1),
                                  latitude = 57.5,
                                  scalefactor = 450.0,
                                  prescribed_velocity = 0.2, 
                                  # pescribed_temperature = t_function,  ####change
                                  # pescribed_salinity = s_function      ####change
                                  )
# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
# CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = s_function, air_concentration = 400)
CO₂_flux = GasExchange(; gas = :CO₂, air_concentration = 400)
# O₂_flux = GasExchange(; gas = :O₂, temperature = t_function, salinity = s_function)
O₂_flux = GasExchange(; gas = :O₂)
model = NonhydrostaticModel(; grid,
                              advection = UpwindBiasedThirdOrder(),  #WENO(;grid)
                              timestepper = :RungeKutta3,
                              closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_photosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          oxygen = true,
                                                          variable_redfield = true,
                                                          open_bottom = true,
                                                          particles = kelp_particles),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     O₂ = FieldBoundaryConditions(top = O₂_flux), ),
                              auxiliary_fields = (; T, S)  ####change
                            #   particles = kelp_particles
                              )

Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002           #in mmolN m^-3
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008           #in mmolN m^-3
sPONᵢ(x, y, z)=0                                                #in mmolN m^-3    Dᵢ
bPONᵢ(x, y, z)=0                                               #in mmolN m^-3     DDᵢ
sPOCᵢ(x, y, z)=0                                                #in mmolC m^-3    Dᶜᵢ
bPOCᵢ(x, y, z)=0                                               #in mmolC m^-3     DDᶜᵢ
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                #in mmolN m^-3
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in mmolN m^-3
DONᵢ(x, y, z)= 0                                             #in mmolN m^-3
DOCᵢ(x, y, z)= 0                                             #in mmolC m^-3
DICᵢ(x, y, z)= 2200                                          #in mmolC m^-3
Alkᵢ(x, y, z)= 2400                                          #in mmolN m^-3
O₂ᵢ(x, y, z) = 16.2*tanh((z+200)/200)+256.6                  #in mmolO m^-3

#set!(model, P=0.03, Z=0.03, NO₃=11.0, NH₄=0.05, DIC=2200.0, Alk=2400.0, O₂=240.0)
set!(model, P=Pᵢ, Z=Zᵢ, sPON=sPONᵢ, bPON=bPONᵢ, sPOC=sPOCᵢ, bPOC=bPOCᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DON=DONᵢ, DOC=DOCᵢ, DIC=DICᵢ, Alk=Alkᵢ, O₂=O₂ᵢ)  #没有速度了
#julia> keys(pp["timeseries"]) 15-element Vector{String}:  "NO₃"  "NH₄" "P" "Z" "sPON" "bPON" "DON" "DIC" "Alk" "O₂" "sPOC" "bPOC" "DOC" "PAR" "t"#

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Couples the particles to the biodeochemical model
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=3.5minutes, stop_time=duration) 

# simulation.callbacks[:couple_particles] = Callback(Particles.infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "kelp450_11_new4"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, model.tracers, filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true) #merge(model.tracers, model.auxiliary_fields),
simulation.output_writers[:particles] = JLD2OutputWriter(model, (; kelp_particles), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)


#simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.1, c_adv=0.5, c_diff=0.5, w = 200/day, relaxation=0.95), TimeStepCallsite())

scale_negative_tracers = ScaleNegativeTracers((:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))   #ScaleNegativeTracers(;model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)) changed (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

# plankton_redfield = model.biogeochemistry.phytoplankton_redfield
# scale_negative_carbon_tracers = ScaleNegativeTracers(tracers = (:P, :Z, :DOC, :sPOC, :bPOC, :DIC), scalefactors = (P = plankton_redfield, Z = plankton_redfield, DOC = 1, sPOC = 1, bPOC = 1, DIC = 1))
# simulation.callbacks[:neg2] = Callback(scale_negative_carbon_tracers; callsite = UpdateStateCallsite())


#wizard = TimeStepWizard(cfl = 0.2, diffusive_cfl = 0.2, max_change = 2.0, min_change = 0.5, cell_diffusion_timescale = column_diffusion_timescale, cell_advection_timescale = column_advection_timescale)
#simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([i*years for i=1:12]), prefix=filename) #prefix="kelp_01"
# ## Run!
# Finally we run the simulation
run!(simulation)


















