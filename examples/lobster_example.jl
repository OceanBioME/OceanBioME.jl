"""
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model and idealized mixed layer depth and surface PAR timeseries
In this case, the script uses an idealized mixed layer depth to construct an idealized diffusivity profile. 
The script uses a timeseries of the photosynthetically available radiation for forcing the LOBSTER model.
This script also shows how to interpolate PAR timeseries data to access to arbitrary time. The PAR profile, as an auxiliary field, will be updated every 1 timestep.
References
(1) Lévy, M., Gavart, M., Mémery, L., Caniaux, G. and Paci, A., 2005. A four‐dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model. Journal of Geophysical Research: Oceans, 110(C7).
(2) Lévy, M., Klein, P. and Treguier, A.M., 2001. Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime. Journal of marine research, 59(4), pp.535-565.
(3) Resplandy, L., Martin, A.P., Le Moigne, F., Martin, P., Aquilina, A., Mémery, L., Lévy, M. and Sanders, R., 2012. How does dynamical spatial variability impact 234Th-derived estimates of organic export?. Deep Sea Research Part I: Oceanographic Research Papers, 68, pp.24-45.
(4) Morel, A. and Maritorena, S., 2001. Bio‐optical properties of oceanic waters: A reappraisal. Journal of Geophysical Research: Oceans, 106(C4), pp.7163-7180.
(5) Resplandy, L., Lévy, M., d'Ovidio, F. and Merlivat, L., 2009. Impact of submesoscale variability in estimating the air‐sea CO2 exchange: Results from a model study of the POMME experiment. Global Biogeochemical Cycles, 23(1).
"""

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
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

using OceanBioME 

params = LOBSTER.defaults  #load parameters in src/parameters/lobster.jl 

# Step 1: import data
#The mixed layer depth is used to construct an idealized diffusivity profile
time_mld = [0.,50,86,105,280,364] # in days, specify piecewiselinear time nodes, must be float 
mld = [250,420,420,20,20,250] # mixed layer depth in meters 
mld_itp = LinearInterpolation((time_mld)days, mld)  #

#The surface PAR timeseries is used to calculate PAR profile
time_par = [0.,20,120,200,330,364]
par = [5,5,90,90,5,5]
PAR_itp = LinearInterpolation((time_par)days, par)
surface_PAR(t) = PAR_itp(mod(t, 364days))  # will be used as part of Option 2 of PAR_field

# Simulation duration    
duration=2days #2years    

# Define the grid
Lx = 1
Ly = 1
Nx = 1
Ny = 1
Nz = 33 # number of points in the vertical direction
Lz = 600 # domain depth            

grid = RectilinearGrid(
                size=(Nx, Ny, Nz), 
                extent=(Lx, Ly, Lz))

#Initialize a PAR field
PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

#Set up the OceanBioME model with the specified biogeochemical model, light, and boundary conditions
#If optional_sets, like carboates or oxygen, are not needed, just set optional_sets=(). The input of topboundaries has to be a NamedTuple.
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field, optional_sets=(), topboundaries=(none=nothing,))

@info "Setup OceanBioME model"
#Create a function of the turbulent vertical diffusivity. This is an idealized functional form, but the depth of mixing is based on an interpolation to the idealized mixed layer depth
κₜ(x, y, z, t) = 1e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4; 

# Now, create a 'model' to run in Oceananignas
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing =  bgc.forcing,
                            boundary_conditions = bgc.boundary_conditions,
                            auxiliary_fields = (PAR=PAR_field, )
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

## set the initial conditions using functions or constants:
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, u=0, v=0, w=0, b=0)

## Set up th simulation
Δt=40 #timestep in second 

simulation = Simulation(model, Δt=Δt, stop_time=duration) 

# create a model 'callback' to update the light (PAR) profile every 1 timestep and integrate sediment model
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))#

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

# Create a message to display every 100 timesteps                                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

#setup dictionary of fields
fields = Dict(zip((["$t" for t in bgc.tracers]..., "PAR"), ([getproperty(model.tracers, t) for t in bgc.tracers]..., [getproperty(model.auxiliary_fields, t) for t in (:PAR, )]...)))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="lobster_example.nc", schedule=TimeInterval(1days), overwrite_existing=true)

@info "Setup simulation"
# Run the simulation                            
run!(simulation)

# Load and plot the results
results = OceanBioME.Plot.load_tracers(simulation)
plot(OceanBioME.Plot.profiles(results)...)

# Save the plot to a PDF file
savefig("lobster_example.pdf")