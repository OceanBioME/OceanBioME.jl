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
using Printf, Oceananigans
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

using OceanBioME 

# Load parameters from src/Models/Biogeochemistry/LOBSTER.jl
params = LOBSTER.defaults  

# Define temperature and salinity as functions of x, y, z, and t(in seconds). The temperature and salinity functions are needed to calculate the air-sea CO2 flux.
t_function(x, y, z, t) = 2.4*cos(t*2π/year + 50day) + 10
s_function(x, y, z, t) = 35.0
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
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

#sediment bcs
sediment=Boundaries.Sediments.Soetaert.setupsediment(grid)

# Set up the OceanBioME model with the specified biogeochemical model, grid, parameters, light, and boundary conditions
bgc = Setup.Oceananigans(:LOBSTER, grid, params, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), bottomboundaries=sediment.boundary_conditions, sinking=true, open_bottom=true)
@info "Setup BGC model"

# Function of the turbulent vertical diffusivity. This is an idealized functional form, and the depth of mixing is based on an analtytical approximation
H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 

model = NonhydrostaticModel(
                                advection = WENO(;grid),
                                timestepper = :RungeKutta3,
                                grid = grid,
                                tracers = bgc.tracers,
                                closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                                forcing =  merge(bgc.forcing, sediment.forcing),
                                boundary_conditions = bgc.boundary_conditions,
                                auxiliary_fields = merge((; PAR), sediment.auxiliary_fields)
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
set!(model, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, DIC=DICᵢ, ALK=ALKᵢ, OXY=OXYᵢ)

## Set up the simulation
simulation = Simulation(model, Δt=5minutes, stop_time=100days)

# create a model 'callback' to update the light (PAR) profile every 1 timestep and integrate sediment model
simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)));

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))
#update the timestep length each day
@warn "Timestep utility may cause instability"
simulation.callbacks[:timestep] = Callback(update_timestep!, TimeInterval(1hour), parameters = (w=200/day, c_diff = 0.45, c_adv = 0.45, relaxation=0.75, c_forcing=1.0, Δt_max=3minute)) 

#setup dictionary of fields
fields = Dict(zip((["$t" for t in bgc.tracers]..., "PAR"), ([getproperty(model.tracers, t) for t in bgc.tracers]..., model.auxiliary_fields.PAR)))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="sediment.nc", schedule=TimeInterval(1days), overwrite_existing=true)

# Oceananians storage of sliced fields is currently broken (https://github.com/CliMA/Oceananigans.jl/issues/2770) so here is a work around
using JLD2

function store_sediment!(sim)
    jldopen("sediment.jld2", "a+") do file
        file["Nᵣ/$(sim.model.clock.time)"] = sim.model.auxiliary_fields.Nᵣ[1, 1, 1]
        file["Nᵣᵣ/$(sim.model.clock.time)"] = sim.model.auxiliary_fields.Nᵣᵣ[1, 1, 1]
        file["Nᵣₑ/$(sim.model.clock.time)"] = sim.model.auxiliary_fields.Nᵣₑ[1, 1, 1]
    end
end

simulation.callbacks[:save_sediment] = Callback(store_sediment!, TimeInterval(1days))

#=
# This is currently broken in Oceananigans 
sediment_fields = Dict(zip(("Nᵣᵣ", "Nᵣ", "Nᵣₑ"), model.auxiliary_fields[(:Nᵣᵣ, :Nᵣ, :Nᵣₑ)]))
simulation.output_writers[:sediment_profiles] = NetCDFOutputWriter(model, sediment_fields, filename="sediment_sediment.nc", schedule=TimeInterval(1days), overwrite_existing=true, indices=(:, :, 1:1))
=#

#doesn't work yet
#simulation.callbacks[:neg_sed] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:Nᵣᵣ, :Nᵣ, :Nᵣₑ), warn=true))

@info "Setup simulation"
ΣN₀ = Budget.calculate_budget(model, true, (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM))
# Run the simulation                            
run!(simulation)

ΣN₁ = Budget.calculate_budget(model, true, (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM))
@info "Nitrogen budget varied from $ΣN₀ to $ΣN₁"
#=
# Load and plot the results
include("PlottingUtilities.jl")
results = load_tracers(simulation)
plot(profiles(results)...)

# Save the plot to a PDF file
savefig("sediment.pdf")
=#