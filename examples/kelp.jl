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

# ## Kelp Particle setup
@info "Setting up kelp particles"
n = 5 # number of kelp fronds
z₀ = [-21:5:-1;]*1.0 # depth of kelp fronds

kelp_particles = SLatissima.setup(;n, 
                                  x₀ = Lx/2, y₀ = Ly/2, z₀, 
                                  A₀ = 30.0, N₀ = 0.01, C₀ = 0.1, 
                                  latitude = 57.5,
                                  scalefactor = 500.0, 
                                  T = t_function, S = s_function, urel = 0.2, 
                                  optional_tracers = (:NH₄, :DIC, :DD, :DDᶜ, :OXY, :DOM))

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))
# ## Biogeochemical and Oceananigans model
# Here we instantiate the LOBSTER model which will return all of the information we then need to pass onto Oceananigans and set the initial conditions.
# > We are being a bit picky about the Oceananigans model setup (e.g. specifying the advection scheme) as this gives the best simple results but you may find differently.
bgc = Setup.Oceananigans(:LOBSTER, grid, params, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc)) 

@info "Setup BGC model"

model = NonhydrostaticModel(
                                                advection = WENO(;grid),
                                                timestepper = :RungeKutta3,
                                                grid = grid,
                                                tracers = bgc.tracers,
                                                closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                                                forcing = bgc.forcing,
                                                boundary_conditions = bgc.boundary_conditions,
                                                auxiliary_fields = (; PAR),
                                                particles = kelp_particles
)
set!(model, P=0.03, Z=0.03, D=0.0, DD=0.0, Dᶜ=0.0, DDᶜ=0.0, NO₃=11, NH₄=0.05, DOM=0.0, DIC=2200.0, ALK=2400.0, OXY=240.0)


# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Couples the particles to the biodeochemical model
# - Update the PAR field from the surface PAR and phytoplankton concentration
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=10minutes, stop_time=100days) 

simulation.callbacks[:couple_particles] = Callback(Particles.infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())

simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(params, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))                
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "kelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))
simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.5, c_adv=0.6, c_diff=0.6, w = 200/day, relaxation=0.75), TimeStepCallsite())

# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visulise the results with some post processing to diagnose the air-sea CO₂ flux - hopefully this looks different to the example without kelp!

P = FieldTimeSeries("$filename.jld2", "P")
NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
Z = FieldTimeSeries("$filename.jld2", "Z")
D = FieldTimeSeries("$filename.jld2", "D") 
DD = FieldTimeSeries("$filename.jld2", "DD")
DIC = FieldTimeSeries("$filename.jld2", "DIC")
Dᶜ = FieldTimeSeries("$filename.jld2", "Dᶜ")
DDᶜ = FieldTimeSeries("$filename.jld2", "DDᶜ")
ALK = FieldTimeSeries("$filename.jld2", "ALK")

x, y, z = nodes(P)
times = P.times

air_sea_CO₂_flux = zeros(size(P)[4])
carbon_export = zeros(size(P)[4])
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = Boundaries.airseaflux(0.0, 0.0, t, DIC[1, 1, Nz, i], ALK[1, 1, Nz, i], t_function(1, 1, 0, t), s_function(1, 1, 0, t),  merge(Boundaries.defaults.airseaflux, (T=t_function, S=s_function, gas=:CO₂)))
    carbon_export[i] = (Dᶜ[1, 1, end-20, i]*LOBSTER.D_sinking .+ DDᶜ[1, 1, end-20, i]*LOBSTER.DD_sinking)
end

using GLMakie
f=Figure(backgroundcolor=RGBf(1, 1, 1), fontsize=30)

axP = Axis(f[1, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Phytoplankton concentration (mmol N/m³)")
hmP = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(P[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbP = Colorbar(f[1, 3], hmP)

axNO₃ = Axis(f[1, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Nitrate concentration (mmol N/m³)")
hmNO₃ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(NO₃[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbNO₃ = Colorbar(f[1, 6], hmNO₃)

axZ = Axis(f[2, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Zooplankton concentration (mmol N/m³)")
hmZ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(Z[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbZ = Colorbar(f[2, 3], hmZ)

axD = Axis(f[2, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Detritus concentration (mmol N/m³)")
hmD = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(D[1, 1, end-23:end, 1:101])' .+ float.(DD[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbD = Colorbar(f[2, 6], hmD)

axfDIC = Axis(f[3, 1:4], xlabel="Time (days)", title="Air-sea CO₂ flux and Sinking", ylabel="Flux (kgCO₂/m²/year)")
hmfDIC = GLMakie.lines!(times./days, air_sea_CO₂_flux.*(12+16*2).*year/(1000*1000), label="Air-sea flux")
hmfExp = GLMakie.lines!(times./days, carbon_export.*(12+16*2).*year/(1000*1000), label="Sinking export")

f[3, 5] = Legend(f, axfDIC, "", framevisible = false)

save("examples/$(filename).png", f)

# We can also have a look at how the kelp particles evolve

include("PlottingUtilities.jl")

res_kelp = load_particles("$(filename)_particles.jld2")

xs, ys = res_kelp.t/(1day), res_kelp.results[3,:,1]

f = Figure(backgroundcolor = RGBf(1, 1, 1), resolution = (1920, 1050))
gInput = f[1, 1] = GridLayout()
gProperties = f[1, 2] = GridLayout()
gOutput = f[1, 3] = GridLayout()

ax1, hm1 = GLMakie.heatmap(gInput[1, 1], xs, ys, res_kelp.results[17,:,:]')
ax1.title = "NO₃"
ax1.xticklabelsvisible= false
cb1 = Colorbar(gInput[1, 1:2], hm1, label = "mmol N/m³")

ax2, hm2 = GLMakie.heatmap(gInput[2, 1], xs, ys, res_kelp.results[18,:,:]')
ax2.title = "NH₄"
ax2.ylabel = "depth (m)"
ax2.xticklabelsvisible= false
cb2 = Colorbar(gInput[2, 1:2], hm2, label = "mmol N/m³")

ax3, hm3 = GLMakie.heatmap(gInput[3, 1], xs, ys, res_kelp.results[19,:,:]')
ax3.title = "PAR"
ax3.xlabel = "time (day)"
cb3 = Colorbar(gInput[3, 1:2], hm3, label = "einstein/m²/day")

ax1, hm1 = GLMakie.heatmap(gProperties[1, 1], xs, ys, res_kelp.results[7,:,:]')
ax1.title = "Area"
ax1.xticklabelsvisible= false
cb1 = Colorbar(gProperties[1, 1:2], hm1, label = "dm²")

ax2, hm2 = GLMakie.heatmap(gProperties[2, 1], xs, ys, ((res_kelp.results[8,:,:].+SLatissima.defaults.N_struct).*SLatissima.defaults.K_A .*res_kelp.results[7, :, :])')
ax2.title = "Total Nitrogen (structural + reserve)"
ax2.ylabel = "depth (m)"
ax2.xticklabelsvisible= false
cb2 = Colorbar(gProperties[2, 1:2], hm2, label = "gN/frond")

ax3, hm3 = GLMakie.heatmap(gProperties[3, 1], xs, ys, ((res_kelp.results[9,:,:].+SLatissima.defaults.C_struct).*SLatissima.defaults.K_A .*res_kelp.results[7, :, :])')
ax3.title = "Total Carbon (structural + reserve)"
ax3.xlabel = "time (day)"
cb3 = Colorbar(gProperties[3, 1:2], hm3, label = "gC/frond")

ax1, hm1 = GLMakie.heatmap(gOutput[1, 1], xs, ys, res_kelp.results[10,:,:]')
ax1.title = "NO₃ uptake"
cb1 = Colorbar(gOutput[1, 1:2], hm1, label = "mmol N/s")
ax1.xticklabelsvisible= false

ax2, hm2 = GLMakie.heatmap(gOutput[2, 1], xs, ys, res_kelp.results[11,:,:]')
ax2.title = "NH₄ uptake"
ax2.xticklabelsvisible= false
cb2 = Colorbar(gOutput[2, 1:2], hm2, label = "mmol N/s")

ax3, hm3 = GLMakie.heatmap(gOutput[3, 1], xs, ys, res_kelp.results[12,:,:]')
ax3.title = "Primary production (photosynthesis - respiration)"
ax3.xticklabelsvisible= false
ax3.ylabel = "depth (m)"
cb3 = Colorbar(gOutput[3, 1:2], hm3, label = "mmol C/s")

ax4, hm4 = GLMakie.heatmap(gOutput[4, 1], xs, ys, res_kelp.results[14,:,:]')
ax4.title = "Exudation (DOC output)"
ax4.xticklabelsvisible= false
cb4 = Colorbar(gOutput[4, 1:2], hm4, label = "mmol C/s")

ax5, hm5 = GLMakie.heatmap(gOutput[5, 1], xs, ys, res_kelp.results[16,:,:]')
ax5.title = "Frond errosion (POM output)"
ax5.xlabel = "time (day)"
cb5 = Colorbar(gOutput[5, 1:2], hm5, label = "mmol C/s")
save("examples/$(filename)_particles.png")
