"""
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model.
In this case, the script 
*uses an idealized mixed layer depth (mld) to construct an idealized diffusivity profile;
*uses an idealized surface photosynthetically available radiation (PAR) timeseries to derive the PAR profile for forcing the LOBSTER model, according to a light absorption model.
    The PAR profile will be updated every 1 timestep as an auxiliary field; and
*shows how to interpolate a timeseries of data to access to arbitrary time. 
References
(1) Lévy, M., Gavart, M., Mémery, L., Caniaux, G. and Paci, A., 2005. A four‐dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model. Journal of Geophysical Research: Oceans, 110(C7).
(2) Lévy, M., Klein, P. and Treguier, A.M., 2001. Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime. Journal of marine research, 59(4), pp.535-565.
(3) Resplandy, L., Martin, A.P., Le Moigne, F., Martin, P., Aquilina, A., Mémery, L., Lévy, M. and Sanders, R., 2012. How does dynamical spatial variability impact 234Th-derived estimates of organic export?. Deep Sea Research Part I: Oceanographic Research Papers, 68, pp.24-45.
(4) Morel, A. and Maritorena, S., 2001. Bio‐optical properties of oceanic waters: A reappraisal. Journal of Geophysical Research: Oceans, 106(C4), pp.7163-7180.
(5) Resplandy, L., Lévy, M., d'Ovidio, F. and Merlivat, L., 2009. Impact of submesoscale variability in estimating the air‐sea CO2 exchange: Results from a model study of the POMME experiment. Global Biogeochemical Cycles, 23(1).
"""

# Load the needed modules
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

# Load parameters from src/Models/Biogeochemistry/LOBSTER.jl
params = LOBSTER.defaults  

# Create an idealized mixed layer depth data, and interpolate it to access to mixed layer depth at arbitrary time
time_mld = [0.,50,86,105,280,364] # in days, piecewiselinear time nodes, must be float 
mld = [250,420,420,20,20,250] # in meters, mixed layer depth 

# Convert unit from days to seconds, then conduct linear interpolation to access to mixed layer depth at arbitrary time
#'days' is a Float64 constant equal to 24hours (defined in Oceananigans). Useful for increasing the clarity of scripts
mld_itp = LinearInterpolation((time_mld)days, mld)  

# Create an idealized surface PAR timeseries
time_par = [0.,20,120,200,330,364] # in days, piecewiselinear time nodes, must be float 
par = [5,5,90,90,5,5] # in W/m^2
PAR_itp = LinearInterpolation((time_par)days, par) #conduct linear interpolation to access to PAR at arbitrary time

# Define surface_PAR as a function of time(in seconds). It will be used to update PAR profile as an auxiliary field
surface_PAR(t) = PAR_itp(mod(t, 364days))  # the remainder of t after floored division by 364days. It creates an annual cycle representation of PAR. 

# Simulation duration    
duration=100days #2years as another example    

# Define the grid
Lx = 20
Ly = 20
Nx = 1
Ny = 1
Nz = 33 # number of points in the vertical direction
Lz = 600 # domain depth            

grid = RectilinearGrid(
                size=(Nx, Ny, Nz), 
                extent=(Lx, Ly, Lz)) #A regular rectilinear grid in this case. 

# Initialize a PAR field
PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# Set up the OceanBioME model with the specified biogeochemical model, grid, parameters, light, and optional boundary conditions(More examples will be provided.)
# If optional tracers are included, add e.g. optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc). 
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field) 

@info "Setup OceanBioME model"
# Create a function of the turbulent vertical diffusivity. This is an idealized functional form, and the depth of mixing is based on an interpolation to the idealized mixed layer depth
κₜ(x, y, z, t) = 1e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4; 

# Create a 'model' to run in Oceananignas
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing = bgc.forcing,
                            boundary_conditions = bgc.boundary_conditions,
                            auxiliary_fields = (PAR=PAR_field, )
                            )

# Initialize the biogeochemical variables
# These initial conditions are set rather arbitrarily in the hope that the model will converge to a repeatable annual cycle if run long enough
Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002           #in mmolN m^-3
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008           #in mmolN m^-3
Dᵢ(x, y, z)=0                                                #in mmolN m^-3
DDᵢ(x, y, z)=0                                               #in mmolN m^-3
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                #in mmolN m^-3
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in mmolN m^-3 
DOMᵢ(x, y, z)= 0                                             #in mmolN m^-3  

# Set the initial conditions using functions or constants:
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, u=0, v=0, w=0, b=0)

# Set up the simulation
Δt=40 #timestep in second 
simulation = Simulation(model, Δt=Δt, stop_time=duration) 

# Create a model 'callback' to update the light (PAR) profile every 1 timestep and integrate sediment model
simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=surface_PAR,))));

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

# Create a message to display every 100 timesteps                                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

# Setup dictionary of fields
fields = Dict(zip((["$t" for t in bgc.tracers]..., "PAR"), ([getproperty(model.tracers, t) for t in bgc.tracers]..., [getproperty(model.auxiliary_fields, t) for t in (:PAR, )]...)))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="column.nc", schedule=TimeInterval(1days), overwrite_existing=true)

@info "Setup simulation"
# Run the simulation                            
run!(simulation)

include("PlottingUtilities.jl")
# Load and plot the results
results = load_tracers(simulation)
plot(profiles(results)...)

# Save the plot to a PDF file
savefig("column.pdf")