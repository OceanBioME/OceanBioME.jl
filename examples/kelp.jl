"""
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model using a PAR timeseries from the North Atlantic subpolar gyre
This script also shows how to read in data files (in netCDF format) for forcing the model
In this case, the script reads in the mixed layer depth from the Mercator Ocean state esimate and uses this to construct an idealized diffusivity profile
The script also reads in a timeseries of the photosynthetically available radiation for forcing the LOBSTER model
References
(1) Lévy, M., Gavart, M., Mémery, L., Caniaux, G. and Paci, A., 2005. A four‐dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model. Journal of Geophysical Research: Oceans, 110(C7).
(2) Lévy, M., Klein, P. and Treguier, A.M., 2001. Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime. Journal of marine research, 59(4), pp.535-565.
(3) Resplandy, L., Martin, A.P., Le Moigne, F., Martin, P., Aquilina, A., Mémery, L., Lévy, M. and Sanders, R., 2012. How does dynamical spatial variability impact 234Th-derived estimates of organic export?. Deep Sea Research Part I: Oceanographic Research Papers, 68, pp.24-45.
(4) Morel, A. and Maritorena, S., 2001. Bio‐optical properties of oceanic waters: A reappraisal. Journal of Geophysical Research: Oceans, 106(C4), pp.7163-7180.
(5) Resplandy, L., Lévy, M., d'Ovidio, F. and Merlivat, L., 2009. Impact of submesoscale variability in estimating the air‐sea CO2 exchange: Results from a model study of the POMME experiment. Global Biogeochemical Cycles, 23(1).
"""
@warn "Not currently known to be working"
# First, load the needed modules
using Random   
using Printf
using Plots
using JLD2
using NetCDF
using HDF5
using Interpolations
using Statistics 

using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

using OceanBioME 

params = LOBSTER.defaults  #load parameters in src/parameters/lobster.jl 

# Step 1: import data
# The temperature and salinity are needed to calculate the air-sea CO2 flux.  The mixed layer depth is used to construct an idealized diffusivity profile
filename = "./OceanBioME_example_data/subpolar.nc" 
time = ncread(filename, "time")
temp = ncread(filename, "temp")
salinity = ncread(filename, "so")
mld = ncread(filename, "mld")
par = ncread(filename, "par")

temp_itp = LinearInterpolation(time, temp) 
salinity_itp = LinearInterpolation(time, salinity) 
mld_itp = LinearInterpolation(time, mld) 
PAR_itp = LinearInterpolation(time, par)

t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))
surface_PAR(t) = PAR_itp(mod(t, 364days))  # will be used as part of Option 2 of PAR_field

# Simulation duration    
duration=2years    

# Define the grid
Lx = 500
Ly = 500
Nx = 1
Ny = 1
Nz = 33 # number of points in the vertical direction
Lz = 600 # domain depth             # subpolar mixed layer depth max 427m 

refinement = 10 # controls spacing near surface (higher means finer spaced)  #1.2 5 -4.9118e9
stretching = 5.754   # controls rate of stretching at bottom
                
## Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
## Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
## Generating function
z_faces(k) = k == 1 ? -(Lz + 4) : Lz * (ζ₀(k-1) * Σ(k-1) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz+1), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)

PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid) #initialize a PAR field 
PAR_func(x, y, z, t) = PAR_extrap(z, mod(t, 364days))    # Define the PAR as a function. z goes first, then t. 

# Specify the boundary condition for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

#Specify the sediment model boundary boundary_conditions
#The sediment model can be turned of by removing the referances to it in the bgc = Setup ... line, and removing the auxillary fields (except PAR) when the Oceananigans model is setup
sediment_bcs=Boundaries.setupsediment(grid)

#Set up the OceanBioME model with the specified biogeochemical model, light, and boundary conditions
@info "Setting up biogeochemical model"
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), bottomboundaries = sediment_bcs.boundary_conditions)

# Setup the kelp particles 
# The first year of model run time won't have any kelp, so we set the area to 0.  We will re-initialize the kelp model after 1 year of model time
@info "Setting up kelp particles"
n_kelp=100 # number of kelp fronds
z₀ = [-100:-1;]*1.0 # depth of kelp fronds
kelp_particles = SLatissima.setup(n_kelp, Lx/2, Ly/2, z₀, 0.0, 0.0, 0.0, 57.5, 100.0, t_function, s_function, 0.2)

# create a function with the vertical turbulent vertical diffusivity. This is an idealized functional form, but the depth of mixing is based on an interpolation to the mixed layer depth from the Mercator Ocean state estimate
κₜ(x, y, z, t) = 1e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation

# Now, create a 'model' to run in Oceananignas
@info "Creating and initalizing Oceananigans model"
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing =  bgc.forcing,
                            boundary_conditions = bgc.boundary_conditions,
                            auxiliary_fields = merge((PAR=PAR_field, ), sediment_bcs.auxiliary_fields),  #comment out this line if using functional form of PAR
                            particles = kelp_particles
                            )

# Initialize the biogeochemical variables
# These initial conditions are set rather arbitrarily in the hope that the model will converge to a repeatable annual cycle if run long enough
Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002         
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008          
Dᵢ(x, y, z)=0
DDᵢ(x, y, z)=0
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4  
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05      
DOMᵢ(x, y, z)= 0 
DICᵢ(x, y, z)= 2200 
ALKᵢ(x, y, z)= 2400
OXYᵢ(x, y, z) = 240

## set the initial conditions using functions or constants:
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, DIC=DICᵢ, ALK=ALKᵢ, OXY=OXYᵢ, u=0, v=0, w=0, b=0)

## Set up the simulation
@info "Setting up the simulation"
#Find the maximum timestep (inefficient as this is only required later on (around 100 days in) when sinking particles speed up)
#Timestep wizard can not be used in this instance as the diffusivity is functional otherwise would be better to use
#Can do about 90s timestep early on but reduces to about 40s later on
#Not sure advective one is relivant
c_diff = 0.18
#c_adv = 0.05
dz²_κ=[findmin(grid.Δzᵃᵃᶜ[i]^2 ./(κₜ.(0.5,0.5,grid.zᵃᵃᶜ[i],[0:364;])))[1] for i in 1:Nz]
Δt=c_diff*findmin(dz²_κ)[1] #min(c_diff*findmin(dz²_κ)[1], -c_adv*findmin(grid.Δzᵃᵃᶜ)[1]/params.V_dd)

simulation = Simulation(model, Δt=Δt, stop_time=min(duration, 1year)) 
# create a model 'callback' to update the light (PAR) profile every 1 timestep and integrate sediment model
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(merge(params, Light.defaults), (surface_PAR=surface_PAR,)))#comment out if using PAR as a function, PAR_func
simulation.callbacks[:integrate_sediment] = sediment_bcs.callback

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

# Create a message to display every 100 timesteps                                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

# Specify which data to save 
#setup dictionary of fields (jld2 saving does not work for irregular grids)
fields = Dict(zip((["$t" for t in bgc.tracers]..., "PAR", "Nᵣ", "Nᵣᵣ", "Nᵣₑ"), ([getproperty(model.tracers, t) for t in bgc.tracers]..., [getproperty(model.auxiliary_fields, t) for t in (:PAR, :Nᵣ, :Nᵣᵣ, :Nᵣₑ)]...)))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="kelp.nc", schedule=TimeInterval(1days))

#checkpoint after warmup so we don't have to rerun for different kelp configs
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([1year]), prefix="kelp_checkpoint")

@info "Running simulation for the first year (without kelp)"
@info "(Note that the first timestep will take some time to complete)"
# Run the simulation                           
run!(simulation)

@info "Initializing kelp"
#reset the kelp properties after warmup
model.particles.properties.A .= 30.0*ones(n_kelp)
model.particles.properties.N .= 0.01*ones(n_kelp)
model.particles.properties.C .= 0.1*ones(n_kelp)

#start recording particles
simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), 
                          filename = "particles.jld2",
                          schedule = TimeInterval(1day),
                          overwrite_existing = true)

simulation.stop_time = duration
simulation.Δt = 65 #seems to be some kind of instability occuring just after kelp is added with the longer timesteps

#run rest of simulation
@info "Restarting simulation to run for the year, now with kelp"
run!(simulation)

# Load and plot the results
results = OceanBioME.Plot.load_tracers(simulation)
profiles = OceanBioME.Plot.profiles(results)
savefig("subpolar.pdf")

particles = OceanBioME.Plot.load_particles(simulation)
particles = OceanBioME.Plot.particles(particles)
savefig("subpolar_kelp_particles.pdf")
