# # Sediment model coupling
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model with a sediment model coupled providing the bottom boundary condition. This shows:
# - How to setup OceanBioME's biogeochemical models
# - How to setup light attenuation
# - How to add a sediment (or other complicated) boundary model
# - How to include optional tracer sets (carbonate chemistry and oxygen)
# - How to visulise results
# > WARNING: This sediment model is currently wrong and will not produce correct results

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
params = LOBSTER.defaults  

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((mod(t, 365days)-200days)/50days)^2))) .+ 2

H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))
MLD(t) = (-10 .-340 .*(1 .-fmld1(364.99999days).*exp.(-t/25days).-fmld1.(mod.(t, 365days))))
κₜ(x, y, z, t) = 1e-2*max(1-(z+MLD(t)/2)^2/(MLD(t)/2)^2,0)+1e-4; 

t_function(x, y, z, t) = 2.4*cos(t*2π/year + 50day) + 10
s_function(x, y, z, t) = 35.0

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 
PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

# For the sediment model we have to setup the tracer field first and then pass these to the sediment model
NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY = CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid)
sediment=Boundaries.Sediments.Soetaert.setup(grid, (;D, DD); POM_w=(D=LOBSTER.D_sinking, DD=LOBSTER.DD_sinking))

# ## Biogeochemical and Oceananigans model
bgc = Setup.Oceananigans(:LOBSTER, grid, params, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), open_bottom=true, bottomboundaries=sediment.boundary_conditions) 

@info "Setup BGC model"

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
set!(model, P=0.03, Z=0.03, D=0.0, DD=0.0, Dᶜ=0.0, DDᶜ=0.0, NO₃=11, NH₄=0.05, DOM=0.0, DIC=2200.0, ALK=2400.0, OXY=240.0)


# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

<<<<<<< HEAD
## Set up the simulation
simulation = Simulation(model, Δt=5minutes, stop_time=50days)
=======
simulation = Simulation(model, Δt=10minutes, stop_time=100days) 
>>>>>>> origin/main

simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "sediment"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))
simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.5, c_adv=0.6, c_diff=0.6, w = 200/day, relaxation=0.75), TimeStepCallsite())

<<<<<<< HEAD
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
=======
# ## Run!
# Finally we run the simulation
>>>>>>> origin/main
run!(simulation)
