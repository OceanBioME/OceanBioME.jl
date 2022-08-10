"""
This script illustrates how to run OceanBioME as a LOBSTER model
References
(1)2005 A four-dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model
(2)2001 Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime
(3)2012 How does dynamical spatial variability impact 234Th-derived estimates of organic export
(4)2001 Bio-optical properties of oceanic waters: A reappraisal
"""

using Random   # load required modules
using Printf
using Plots
using JLD2
using NetCDF
using HDF5
using Interpolations
using Statistics 

using Oceananigans#9e8cae18-63c1-5223-a75c-80ca9d6e9a09
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

using BGC 

params = LOBSTER.default  #load parameters in src/parameters/lobster.jl 

######## 1: import annual cycle physics data  #Global Ocean 1/12° Physics Analysis and Forecast updated Daily
filename1 = "subpolar_physics.nc" # One may need to modify path if necessary.  
#ncinfo(filename1)          # ncinfo can be used to check out info in .nc files. 
time_series_second = (0:364)days # Start from zero if we don't use extrapolation. We cannot use extrapolation if we want to use annual cycle.  
so = ncread(filename1, "so");  #salinity Practical Salinity Unit
so_scale_factor = ncgetatt(filename1, "so", "scale_factor") #Real_Value = (Display_Value X scale_factor) + add_offset
so_add_offset = ncgetatt(filename1, "so", "add_offset")
salinity = mean(so, dims=(1,2))[1:365]*so_scale_factor.+so_add_offset # use [1:365] because sometimes one year has 366 days. 
salinity_itp = LinearInterpolation(time_series_second, salinity) #converted to interpolations to access them at arbitary time, how to use it: salinity_itp(mod(timeinseconds,364days))
#plot(salinity)
thetao = ncread(filename1, "thetao");  #temperature,  Degrees Celsius
thetao_scale_factor = ncgetatt(filename1, "thetao", "scale_factor") 
thetao_add_offset = ncgetatt(filename1, "thetao", "add_offset")
temperature = mean(thetao, dims=(1,2))[1:365]*thetao_scale_factor.+thetao_add_offset
temperature_itp = LinearInterpolation(time_series_second, temperature)  
#plot(temperature)
mlotst = ncread(filename1, "mlotst"); #mixed_layer_depth,Meters
mlotst_scale_factor = ncgetatt(filename1, "mlotst", "scale_factor") 
mlotst_add_offset = ncgetatt(filename1, "mlotst", "add_offset")
mixed_layer_depth = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset
mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth)  
#plot(mixed_layer_depth)

######## 2: import annual cycle chl data (to calculate PAR_func) #Global Ocean Biogeochemistry Analysis and Forecast
filename2 = "subpolar_chl.nc"    #subpolar_chl.nc
chl = ncread(filename2, "chl");  #chl scale_factor=1, add_offset=0
chl_mean = mean(chl, dims=(1,2))[1,1,:,1:365] # mg m-3, unit no need to change. 
depth_chl = ncread(filename2, "depth");
#heatmap(1:365, -depth_chl[end:-1:1], chl_mean[end:-1:1,:])

######## 3: import annual cycle surface photosynthetic available radiation (PAR) data #Ocean Color  VIIRS-SNPP PAR daily 9km
path="./subpolar/"    #subtropical   #./subpolar/
par_mean_timeseries=zeros(365)
for i in 1:365    #https://discourse.julialang.org/t/leading-zeros/30450
    string_i = lpad(string(i), 3, '0')
    filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc"
    fid = h5open(filename3, "r")
    par=read(fid["level-3_binned_data/par"]) # read PAR data from file
    BinList=read(fid["level-3_binned_data/BinList"])  # read information on PAR bins.   (:bin_num, :nobs, :nscenes, :weights, :time_rec) 
    par_mean_timeseries[i] = mean([par[i][1]/BinList[i][4] for i in 1:length(par)])*3.99e-10*545e12/(1day)  #average PAR values in bins and convert from einstin/m^2/day to W/m^2
end

surface_PAR_itp = LinearInterpolation((0:364)day, par_mean_timeseries)
surface_PAR(t) = surface_PAR_itp(mod(t, 364days))  # will be used as part of Option 2 of PAR_field

#=
Option 1: use PAR as a function, see below PAR_func(x, y, z, t) = PAR_extrap(z, mod(t, 364days))
Obtain PAR_extrap(z, mod(t, 364days)) with annual cycle data chl_mean and surface PAR data par_mean_timeseries as input 
Take the following configs:
* comment out auxiliary_fields = (PAR=PAR_field, ) in Model instantiation model = NonhydrostaticModel(..., #auxiliary_fields = (PAR=PAR_field, ))
* comment out simulation.callbacks[:update_par]
* indicate PAR_func as an input in bgc = Setup.Oceananigans()

Option 2: calculate PAR as a field
Initialize PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)
Update PAR_field through simulation.callbacks[:update_par] and src/Light.jl. Instead of using annual cycle data chl_mean as input, one uses chl as
a function of phytoplankton to calculate PAR. See (A22) in reference (2).
Take the following configs:
* keep auxiliary_fields = (PAR=PAR_field, ) in Model instantiation model = NonhydrostaticModel(..., #auxiliary_fields = (PAR=PAR_field, ))
* keep simulation.callbacks[:update_par]
* indicate PAR_field as an input in bgc = Setup.Oceananigans()
=#

PAR = zeros(length(depth_chl),365)   # initialize PAR
PAR_r = zeros(length(depth_chl),365) # break it down into two wavebands: red and blue.
PAR_b = zeros(length(depth_chl),365)
#PAR_g = zeros(length(depth_chl),365)
PAR[1,:] = par_mean_timeseries       # surface PAR was obtained from Ocean Color shown above
PAR_r[1,:] = par_mean_timeseries/2   # visible light is split equally between two bands. 
PAR_b[1,:] = par_mean_timeseries/2
#PAR_g[1,:] = par_mean_timeseries/3
for i =2:length(depth_chl)  # equations refer to the APPENDIX of reference (2) above 
    PAR_r[i,:]=PAR_r[i-1,:].*exp.(-(params.k_r0.+params.Χ_rp*(chl_mean[i-1,:]).^params.e_r)*(depth_chl[i]-depth_chl[i-1]))
    PAR_b[i,:]=PAR_b[i-1,:].*exp.(-(params.k_b0.+params.Χ_bp*(chl_mean[i-1,:]).^params.e_b)*(depth_chl[i]-depth_chl[i-1]))
    #PAR_g[i,:]=PAR_g[i-1,:].*exp.(-(params.k_g0.+params.Χ_gp*(chl_mean[i-1,:]).^params.e_g)*(depth_chl[i]-depth_chl[i-1]))
    #PAR[i,:]=params.β₁*PAR_b[i,:]+params.β₂*PAR_g[i,:]+params.β₃*PAR_r[i,:]
    PAR[i,:]=PAR_b[i,:]+PAR_r[i,:]
end
#heatmap(1:365, -depth_chl[end:-1:1], PAR[end:-1:1,:])
PAR_itp = Interpolations.interpolate((-depth_chl[end:-1:1], (0:364)day), PAR[end:-1:1,:], Gridded(Linear())) #interpolate to access to them at arbitary time and depth
PAR_extrap = extrapolate(PAR_itp, (Line(),Throw()))  #  PAR_extrap(z, mod(t,364days))  Interpolations.extrapolate Method

# Simulation duration    
duration=1days    #2years

# Define the grid
Lx = 1   #500
Ly = 500
Nx = 1
Ny = 1
Nz = 150#150 # number of points in the vertical direction
Lz = 600 # domain depth             # subpolar mixed layer depth max 427m 


refinement = 5 # controls spacing near surface (higher means finer spaced)  #1.2 5 -4.9118e9
stretching = 2.5   # controls rate of stretching at bottom    #24 2.5  5.754
                
## Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz
## Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement
## Bottom-intensified stretching function 
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
## Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)
#=
grid = RectilinearGrid(
                size=(Nx, Ny, Nz), 
                extent=(Lx, Ly, Lz))
=#

B₀ = 0e-8    #m²s⁻³     #buoyancy tracer 
N² = 9e-6    #dbdz=N^2, s⁻²
buoyancy_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(B₀),
                                       bottom = GradientBoundaryCondition(N²))

########## u boundary condition                                
u₁₀ = 0     # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3  # dimensionless drag coefficient
ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
ρₐ = 1.225   # kg m⁻³, average density of air at sea-level
Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) .+ 273.15  # will be used to calculate air-sea-flux 
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))               # will be used to calculate air-sea-flux 

PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid) #initialize a PAR field 
PAR_func(x, y, z, t) = PAR_extrap(z, mod(t, 364days))    #LoadError: cannot define function PAR; it already has a value. So I renamed it as PAR_func. z goes first, then t. 

dic_bc = Boundaries.setupdicflux(params; forcings=(T=t_function, S=s_function)) #calculate air-sea-flux using DIC, ALK, temperature and salinity. 
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR_field, topboundaries=(DIC=dic_bc, ), optional_sets=(:carbonates, )) #input either PAR_func or PAR_field

#npz = Setup.Oceananigans(:NPZ, grid, NPZ.defaults)
#κₜ(x, y, z, t) = 1e-2*max(1-(z+50)^2/50^2,0)+1e-5;
κₜ(x, y, z, t) = 8e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4; #setup viscosity and diffusivity in the following Model instantiation

###Model instantiation
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing =  bgc.forcing,
                            boundary_conditions = merge((u=u_bcs, b=buoyancy_bcs), bgc.boundary_conditions),
                            auxiliary_fields = (PAR=PAR_field, )  #comment out this line if using functional form of PAR
                            )

## Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise


#set initial conditions
initial_mixed_layer_depth = -100 # m
stratification(z) = z < initial_mixed_layer_depth ? N² * z : N² * (initial_mixed_layer_depth)
bᵢ(x, y, z) = stratification(z)         #

#could initiate P from chl data?
Pᵢ(x,y,z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002          # ((tanh((z+100)/50)-1)/2*0.23+0.23)*16/106  
Zᵢ(x,y,z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008          # ((tanh((z+100)/50)-1)/2*0.3+0.3)*16/106         
Dᵢ(x,y,z)=0
DDᵢ(x,y,z)=0
NO₃ᵢ(x,y,z)= (1-tanh((z+300)/150))/2*6+11.4   #  # 17.5*(1-tanh((z+100)/10))/2
NH₄ᵢ(x,y,z)= (1-tanh((z+300)/150))/2*0.05+0.05       #1e-1*(1-tanh((z+100)/10))/2
DOMᵢ(x,y,z)= 0 
DICᵢ(x,y,z)= 2200   #  mmol/m^-3
ALKᵢ(x,y,z)= 2400   #  mmol/m^-3

## `set!` the `model` fields using functions or constants:
set!(model, b=bᵢ, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ,DIC=DICᵢ,ALK=ALKᵢ, u=0, v=0, w=0)

## Setting up a simulation

simulation = Simulation(model, Δt=50, stop_time=duration)  #Δt=0.5*(Lz/Nz)^2/1e-2,
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))#comment out if using PAR as a function, PAR_func

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

#=
#this can move into LOBSTER.jl at some point
pco2_bc = zeros(2,round(Int,duration/1day)+3)   #Int(duration/simulation.Δt)
function pco2(sim)  #https://clima.github.io/OceananigansDocumentation/stable/generated/baroclinic_adjustment/
    #i+=1
    pco2_bc[2,round(Int,sim.model.clock.time/1day)+1] = LOBSTER.air_sea_flux(1, 1, sim.model.clock.time, model.tracers.DIC[1,1,end-3], model.tracers.ALK[1,1,end-3], params)*(1years/1000)/(7.7e-4*params.U_10^2)+params.pCO2_air   #1×1×150 Field
    pco2_bc[1,round(Int,sim.model.clock.time/1day)+1] = sim.model.clock.time/1day
    #sim.model.clock.iteration
end 

simulation.callbacks[:pco2] = Callback(pco2, IterationInterval(Int(1day/simulation.Δt))) #callback every 1 day 
=#

# We then set up the simulation:

# Vertical slice
simulation.output_writers[:profiles] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, model.auxiliary_fields),
                          filename = "profile_subpolar_test_PARfield.jld2",
                          indices = (1, 1, :),
                          schedule = TimeInterval(1days),     #TimeInterval(1days),
                            overwrite_existing = true)

#simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), 
#                            prefix = "particles",
#                          schedule = TimeInterval(1minute),
#                             force = true)

# We're ready:

run!(simulation)
#jldsave("pco2_water_subpolar.jld2"; pco2_bc)  
#jldopen("pco2_water_subpolar.jld2", "r")


# ## Turbulence visualization
#
# We animate the data saved in `ocean_wind_mixing_and_convection.jld2`.
# We prepare for animating the flow by creating coordinate arrays,
# opening the file, building a vector of the iterations that we saved
# data at, and defining functions for computing colorbar limits:

## Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)
xb, yb, zb = nodes(model.tracers.b)

results = BGC.Plot.load_tracers(simulation)
BGC.Plot.profiles(results)
savefig("subpolar_test_PARfield.pdf")
