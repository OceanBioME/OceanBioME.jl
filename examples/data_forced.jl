"""
This script illustrates how to run OceanBioME as a 1D column model using the LOBSTER biogeochemical model.
In this case, the script 
*introduces two additonal tracer types: carbonates (DIC and ALK) and oxygen (OXY);
*reads in the temperature, salinity, and mixed layer depth from the Mercator Ocean state esimate of the North Atlantic subpolar gyre;
    The temperature and salinity are needed to calculate the air-sea CO2 flux.  
    The mixed layer depth is used to construct an idealized diffusivity profile.
*reads in the surface PAR timeseries data of the same geographic region as above from the NASA's OceanColor Web; 
    Both the data downloaded are saved seperately and a small sample of them are stored in subpolar.nc for ease of use.
    The surface PAR timeseries is used to derive the PAR profile for forcing the LOBSTER model based on a light absorption model.
    The PAR profile will be updated every 1 timestep as an auxiliary field.
*shows how to read *.nc files.
References
(1) Lévy, M., Gavart, M., Mémery, L., Caniaux, G. and Paci, A., 2005. A four‐dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model. Journal of Geophysical Research: Oceans, 110(C7).
(2) Lévy, M., Klein, P. and Treguier, A.M., 2001. Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime. Journal of marine research, 59(4), pp.535-565.
(3) Resplandy, L., Martin, A.P., Le Moigne, F., Martin, P., Aquilina, A., Mémery, L., Lévy, M. and Sanders, R., 2012. How does dynamical spatial variability impact 234Th-derived estimates of organic export?. Deep Sea Research Part I: Oceanographic Research Papers, 68, pp.24-45.
(4) Morel, A. and Maritorena, S., 2001. Bio‐optical properties of oceanic waters: A reappraisal. Journal of Geophysical Research: Oceans, 106(C4), pp.7163-7180.
(5) Resplandy, L., Lévy, M., d'Ovidio, F. and Merlivat, L., 2009. Impact of submesoscale variability in estimating the air‐sea CO2 exchange: Results from a model study of the POMME experiment. Global Biogeochemical Cycles, 23(1).
"""

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
using Oceananigans.Operators: ∂zᶜᶜᶜ
using OceanBioME 

# Load parameters from src/Models/Biogeochemistry/LOBSTER.jl
params = LOBSTER.defaults  

# Import data
# The temperature and salinity are needed to calculate the air-sea CO2 flux.  The mixed layer depth is used to construct an idealized diffusivity profile.
filename = "./OceanBioME_example_data/subpolar.nc" #A small sample of data downloaded is stored in subpolar.nc for ease of use.
time = ncread(filename, "time")    # time in seconds
temp = ncread(filename, "temp")    # temperature in Degrees Celsius 
salinity = ncread(filename, "so")  # salinity in Practical Salinity Unit
mld = ncread(filename, "mld")      # mixed layer depth in Meters
par = ncread(filename, "par")      # photosynthetically available radiation in W/m^2

# Linear interpolation to access temperature, salinity, mld, and surface PAR at arbitrary time
temperature_itp = LinearInterpolation(time, temp) 
salinity_itp = LinearInterpolation(time, salinity) 
mld_itp = LinearInterpolation(time, mld) 
PAR_itp = LinearInterpolation(time, par)

# Define temperature and salinity as functions of x, y, z, and t(in seconds). The temperature and salinity functions are needed to calculate the air-sea CO2 flux.
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) # the remainder of t after floored division by 364days. It creates an annual cycle representation of temperature. 
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))
# Define surface_PAR as a function of time(in seconds). 
surface_PAR(t) = PAR_itp(mod(t, 364days))  # the remainder of t after floored division by 364days. It creates an annual cycle representation of PAR.

# Simulation duration    
duration=1year

# Define the grid
Lx = 20
Ly = 20
Nx = 1
Ny = 1
Nz = 45 # number of points in the vertical direction
Lz = 600 # domain depth

# Generate vertically stretched grid 
refinement = 8 # controls spacing near surface (higher means finer spaced)  
stretching = 5.754   # controls rate of stretching at bottom      
# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
# Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)     

# Initialize a PAR field                       
PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

# Set up the OceanBioME model with the specified biogeochemical model, grid, parameters, light, and boundary conditions
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), open_bottom=false)
@info "Set up OceanBioME model"

# Create a function with the turbulent vertical diffusivity. This is an idealized functional form, but the depth of mixing is based on an interpolation to the mixed layer depth from the Mercator Ocean state estimate.
κₜ(x, y, z, t) = 2e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4;   #m^2 s^-1

# Create a 'model' to run in Oceananignas
model = NonhydrostaticModel(
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
Pᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002           #in mmolN m^-3
Zᵢ(x, y, z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008           #in mmolN m^-3
Dᵢ(x, y, z)=0                                                #in mmolN m^-3
DDᵢ(x, y, z)=0                                               #in mmolN m^-3
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                #in mmolN m^-3
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in mmolN m^-3
DOMᵢ(x, y, z)= 0                                             #in mmolN m^-3
DICᵢ(x, y, z)= 2200                                          #in mmolC m^-3
ALKᵢ(x, y, z)= 2400                                          #in mmolN m^-3
OXYᵢ(x, y, z) = 240                                          #in mmolO m^-3

# Set the initial conditions using functions or constants:
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, DIC=DICᵢ, ALK=ALKᵢ, OXY=OXYᵢ, u=0, v=0, w=0, b=0)

# Set up the simulation
Δt=0.5minutes
simulation = Simulation(model, Δt=Δt, stop_time=duration) 

# Create a model 'callback' to update the light (PAR) profile every 1 timestep
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(merge(params, Light.defaults), (surface_PAR=surface_PAR,)))

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
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="data_forced.nc", schedule=TimeInterval(1days), overwrite_existing=true)

budget = []

function disp_budget(sim)
    push!(budget, sum(sim.model.tracers.NO₃[1, 1, 1:Nz])+sum(sim.model.tracers.NH₄[1, 1, 1:Nz])+sum(sim.model.tracers.P[1, 1, 1:Nz])+sum(sim.model.tracers.Z[1, 1, 1:Nz])+sum(sim.model.tracers.DOM[1, 1, 1:Nz])+sum(sim.model.tracers.D[1, 1, 1:Nz])+sum(sim.model.tracers.DD[1, 1, 1:Nz]))
end
simulation.callbacks[:budget] = Callback(disp_budget, IterationInterval(1))
@info "Setup simulation"

simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (w=200/day, c_diff = 0.45, c_adv = 0.45, relaxation=0.75))

run!(simulation)
#=
include("PlottingUtilities.jl")
# Load and plot the results
results = load_tracers(simulation)
plot(profiles(results)...)

# Save the plot to a PDF file
savefig("data_forced.pdf")
=#
