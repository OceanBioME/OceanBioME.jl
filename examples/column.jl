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
using Printf

using Oceananigans
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

using OceanBioME 

# Load parameters from src/Models/Biogeochemistry/LOBSTER.jl
params = LOBSTER.defaults  

# Surface PAR function
PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((t-200days)/50days)^2))) .+ 2

# Define the grid
Lx = 20
Ly = 20
Nx = 1
Ny = 1
Nz = 50 # number of points in the vertical direction
Lz = 200 # domain depth            

grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz)) #A regular rectilinear grid in this case. 

# Initialize a PAR field
PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# Set up the OceanBioME model with the specified biogeochemical model, grid, parameters, and optionally boundary conditions and other tracer sets (e.g. carbonate chemistry)
bgc = Setup.Oceananigans(:LOBSTER, grid, params; sinking=true, open_bottom=true) 

@info "Setup OceanBioME model"

# Function of the turbulent vertical diffusivity. This is an idealized functional form, and the depth of mixing is based on an analtytical approximation
H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 

# Create a 'model' to run in Oceananignas
model = NonhydrostaticModel(
                                                advection = WENO(;grid),
                                                timestepper = :RungeKutta3,
                                                grid = grid,
                                                tracers = bgc.tracers,
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
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ)

# Set up the simulation
Δt = 10minutes
simulation = Simulation(model, Δt=Δt, stop_time=100days) 

# Create a model 'callback' to update the light (PAR) profile every 1 timestep and integrate sediment model
simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)));

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

# prevents tracers going negative, numerically questionable but positivity preserving time stepper to be implimented soon, see discussion
simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), ))
@info "Setup simulation"
               
run!(simulation)


include("PlottingUtilities.jl")
# Load and plot the results
results = load_tracers("column.nc", (LOBSTER.tracers..., :PAR), 1)
Plots.plot(profiles(results)...)
savefig("column.png")