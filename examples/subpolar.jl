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

# First, load the needed modules
using Random 
using Printf
using Plots
using JLD2
using NetCDF
using HDF5
using Interpolations
using Statistics 
using Adapt 

using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

using OceanBioME 

const params = LOBSTER.defaults  #load parameters in src/parameters/lobster.jl 

# Step 1: import data from the Mercator Ocean state estimate: Global Ocean 1/12° Physics Analysis and Forecast updated Daily
# The temperature and salinity are needed to calculate the air-sea CO2 flux.  The mixed layer depth is used to construct an idealized diffusivity profile
filename1 = "OceanBioME_example_data/subpolar_physics.nc" # Specify the file containing the Mercator model output. One may need to modify path if necessary.  
time_series_second = (0:364)days # Create an array of times, one per day in units of seconds. Start from zero if we don't use extrapolation. We cannot use extrapolation if we want to use annual cycle.  
so = ncread(filename1, "so");  # Read the salinity in Practical Salinity Units
so_scale_factor = ncgetatt(filename1, "so", "scale_factor") #Real_Value = (Display_Value X scale_factor) + add_offset
so_add_offset = ncgetatt(filename1, "so", "add_offset")
const salinity = mean(so, dims=(1,2))[1:365]*so_scale_factor.+so_add_offset # NEEDS EXPLANATION
const salinity_itp = LinearInterpolation(time_series_second, salinity) # Create a function to interpolate the salinity to a specified time.  To use this function, call: salinity_itp(mod(timeinseconds,364days))
thetao = ncread(filename1, "thetao");  # Read the temperature in units of Degrees Celsius
thetao_scale_factor = ncgetatt(filename1, "thetao", "scale_factor") 
thetao_add_offset = ncgetatt(filename1, "thetao", "add_offset")
const temperature = mean(thetao, dims=(1,2))[1:365]*thetao_scale_factor.+thetao_add_offset # NEEDS EXPLANATION
const temperature_itp = LinearInterpolation(time_series_second, temperature) # Create a function to interpolate the temperature to a specified time
mlotst = ncread(filename1, "mlotst"); # Read the mixed layer depth in units of eters
mlotst_scale_factor = ncgetatt(filename1, "mlotst", "scale_factor") 
mlotst_add_offset = ncgetatt(filename1, "mlotst", "add_offset")
const mixed_layer_depth = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset # NEEDS EXPLANATION
const mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth) # Create a function to interpolate the mixed layer depth to a specified time  

# Step 2: import annual cycle chl data (to calculate PAR_func) #Global Ocean Biogeochemistry Analysis and Forecast
filename2 = "OceanBioME_example_data/subpolar_chl.nc"    #subpolar_chl.nc
#chl = ncread(filename2, "chl");  #chl scale_factor=1, add_offset=0
#chl_mean = mean(chl, dims=(1,2))[1,1,:,1:365] # mg m-3, unit no need to change. 
#depth_chl = ncread(filename2, "depth");

# Step 3: import annual cycle surface photosynthetic available radiation (PAR) data #Ocean Color  VIIRS-SNPP PAR daily 9km
path="./OceanBioME_example_data/subpolar/"    #subtropical   #./subpolar/
par_mean_timeseries=zeros(365)
for i in 1:365    #https://discourse.julialang.org/t/leading-zeros/30450
    string_i = lpad(string(i), 3, '0')
    filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc"
    fid = h5open(filename3, "r")
    par=read(fid["level-3_binned_data/par"]) # read PAR data from file
    BinList=read(fid["level-3_binned_data/BinList"])  # read information on PAR bins.   (:bin_num, :nobs, :nscenes, :weights, :time_rec) 
    par_mean_timeseries[i] = mean([par[i][1]/BinList[i][4] for i in 1:length(par)])*3.99e-10*545e12/(1day)  #average PAR values in bins and convert from einstin/m^2/day to W/m^2
end

const surface_PAR_itp = LinearInterpolation((0:364)day, par_mean_timeseries)
surface_PAR(t) = surface_PAR_itp(mod(t, 364days))  # will be used as part of Option 2 of PAR_field

#=
The light (PAR) can be specified in 1 of 2 ways:
Option 1: use PAR as a function, see below PAR_func(x, y, z, t) = PAR_extrap(z, mod(t, 364days))
Obtain PAR_extrap(z, mod(t, 364days)) with annual cycle data chl_mean and surface PAR data par_mean_timeseries as input 
Use the following steps:
* comment out auxiliary_fields = (PAR=PAR_field, ) in Model instantiation model = NonhydrostaticModel(..., #auxiliary_fields = (PAR=PAR_field, ))
* comment out simulation.callbacks[:update_par]
* indicate PAR_func as an input in bgc = Setup.Oceananigans()
Option 2: calculate PAR as a field
Initialize PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)
Update PAR_field through simulation.callbacks[:update_par] and src/Light.jl. Instead of using annual cycle data chl_mean as input, one uses chl as
a function of phytoplankton to calculate PAR. See (A22) in reference (2).
Use the following steps:
* keep auxiliary_fields = (PAR=PAR_field, ) in Model instantiation model = NonhydrostaticModel(..., #auxiliary_fields = (PAR=PAR_field, ))
* keep simulation.callbacks[:update_par]
* indicate PAR_field as an input in bgc = Setup.Oceananigans()
=#
#=
PAR = zeros(length(depth_chl),365)   # initialize PAR
PAR_r = zeros(length(depth_chl),365) # break it down into two wavebands: red and blue.
PAR_b = zeros(length(depth_chl),365)
PAR[1,:] = par_mean_timeseries       # surface PAR was obtained from Ocean Color shown above
PAR_r[1,:] = par_mean_timeseries/2   # visible light is split equally between two bands. 
PAR_b[1,:] = par_mean_timeseries/2
for i =2:length(depth_chl)  # equations refer to the APPENDIX of reference (2) above 
    PAR_r[i,:]=PAR_r[i-1,:].*exp.(-(params.k_r0.+params.Χ_rp*(chl_mean[i-1,:]).^params.e_r)*(depth_chl[i]-depth_chl[i-1]))
    PAR_b[i,:]=PAR_b[i-1,:].*exp.(-(params.k_b0.+params.Χ_bp*(chl_mean[i-1,:]).^params.e_b)*(depth_chl[i]-depth_chl[i-1]))
    PAR[i,:]=PAR_b[i,:]+PAR_r[i,:]
end
PAR_itp = Interpolations.interpolate((-depth_chl[end:-1:1], (0:364)day), PAR[end:-1:1,:], Gridded(Linear())) # create an interpolation function to allow access to the PAR at arbitary time and depth
PAR_extrap = extrapolate(PAR_itp, (Line(),Throw()))  #  PAR_extrap(z, mod(t,364days))  Interpolations.extrapolate Method
=#
# Simulation duration    
const duration=2years    

# Define the grid
Lx = 500
Ly = 500
Nx = 3
Ny = 3
Nz = 33#150 # number of points in the vertical direction
Lz = 600 # domain depth             # subpolar mixed layer depth max 427m 

refinement = 5 # controls spacing near surface (higher means finer spaced)  #1.2 5 -4.9118e9
stretching = 5.754   # controls rate of stretching at bottom    #24 2.5  5.754
                
## Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
## Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
## Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny, Nz), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)

t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))

PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid) #initialize a PAR field 
PAR_func(x, y, z, t) = PAR_extrap(z, mod(t, 364days))    # Define the PAR as a function. z goes first, then t. 

# Specify the boundary condition for DIC and OXY based on the air-sea CO₂ and O₂ flux
const dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
const oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

#Specify the sediment model boundary boundary_conditions
#The sediment model can be turned of by removing the referances to it in the bgc = Setup ... line, and removing the auxillary fields (except PAR) when the Oceananigans model is setup
const sediment_bcs=Boundaries.setupsediment(grid)

#Set up the OceanBioME model with the specified biogeochemical model, light, and boundary conditions
const bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), bottomboundaries = sediment_bcs.boundary_conditions)
@info "Setup BGC model"

# create a function with the vertical turbulent vertical diffusivity. This is an idealized functional form, but the depth of mixing is based on an interpolation to the mixed layer depth from the Mercator Ocean state estimate
κₜ(x, y, z, t) = 8e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation

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
                            auxiliary_fields = merge((PAR=PAR_field, ), sediment_bcs.auxiliary_fields)  #comment out this line if using functional form of PAR
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

## Set up th simulation
#Find the maximum timestep (inefficient as this is only required later on (around 100 days in) when sinking particles speed up)
#Timestep wizard can not be used in this instance as the diffusivity is functional otherwise would be better to use
#91s is about the maximum time step viable with this diffusivity and sinking specified
#The limit is initially diffusive (needing around 40s timestep), but then becomes advective when sinking particles reach the deeper (larger) cells
c_diff = 6.803014480377265e-7
c_adv = 0.00274
dz²_κ=[findmax(adapt(Array, grid.Δzᵃᵃᶜ)[i]^2 ./(κₜ.(0.5,0.5,adapt(Array, grid.zᵃᵃᶜ)[i],[0:364;])))[1] for i in 1:Nz]
Δt=max(c_diff*findmax(dz²_κ)[1], -c_adv*findmax(adapt(Array, grid.Δzᵃᵃᶜ))[1]/params.V_dd)

simulation = Simulation(model, Δt=Δt, stop_time=duration) 

# create a model 'callback' to update the light (PAR) profile every 1 timestep
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))#comment out if using PAR as a function, PAR_func

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

# Create a message to display every 100 timesteps                                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

# Specify which data to save
#= jld2 version does not work for non regular grid spacing
simulation.output_writers[:profiles] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, model.auxiliary_fields),
                          filename = "profile.jld2",
                          indices = (1, 1, :),
                          schedule = TimeInterval(1days),     #TimeInterval(1days),
                            overwrite_existing = true)=#

#setup dictionary of fields
#fields = Dict(zip((["$t" for t in bgc.tracers]..., "PAR", "Nᵣ", "Nᵣᵣ", "Nᵣₑ"), ([getproperty(model.tracers, t) for t in bgc.tracers]..., [getproperty(model.auxiliary_fields, t) for t in (:PAR, :Nᵣ, :Nᵣᵣ, :Nᵣₑ)]...)))
#simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="subpolar.nc", schedule=TimeInterval(1days), overwrite_existing=true)

@info "Setup simulation"
# Run the simulation                            
run!(simulation)

# Load and plot the results
results = OceanBioME.Plot.load_tracers(simulation)
OceanBioME.Plot.profiles(results)

# Save the plot to a PDF file
savefig("subpolar.pdf")
