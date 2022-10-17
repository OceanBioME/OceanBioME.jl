using OceanBioME
using Oceananigans
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Fields: Field
using Oceananigans.TurbulenceClosures: DiscreteDiffusionFunction
using Printf

#const ϕ = 55
ϕ = 55
@inline function L_day(t)
    #see https://www.sciencedirect.com/science/article/pii/030438009400034F for more complicated options
    θ = 23.45*sin(2π*(283+mod(t/day, 365))/365)*π/180
    ζ = -tan(ϕ*π/180)*tan(θ)
    ζ = ifelse(abs(ζ) > 1, sign(ζ), ζ)
    hA = acos(ζ)*180/π
    return hA/(12*15)
end

params = merge(PISCES.defaults, (; ϕ, L_day)) 

#= Define the grid
Lx = 100
Ly = 100
Lz = 1000 timestep
Nx = 5
Ny = 5
Nz = 40 

# Generate vertically stretched grid 
refinement = 10 # controls spacing near surface (higher means finer spaced)  
stretching = 5.754   # controls rate of stretching at bottom      
h(k) = (k - 1) / Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(size = (Nx, Ny, Nz), 
                                    x = (0, Lx),
                                    y = (0, Ly),
                                    z = z_faces)=#

grid = RectilinearGrid(size=(1, 1, 40), extent=(10, 10, 200))

# Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
#dic_bc = Boundaries.airseasetup(:CO₂)
#oxy_bc = Boundaries.airseasetup(:O₂)

# TODO: Add sediment model (of either PISCES formulation, or option to relax back to external model values)

# Set up the OceanBioME model with the specified biogeochemical model, grid, parameters, light, and boundary conditions
bgc = Setup.Oceananigans(:PISCES, grid, params)#, topboundaries=(DIC=dic_bc, OXY=oxy_bc))

κ(i, j, k, grid, clock, model_fields, params) = params.κ₀#*max(1 - (grid.zᵃᵃᶜ[k] - model_fields.zₘₓₗ[i, j])/2)^2 / (-model_fields.zₘₓₗ[i, j])/2)^2,0)+1e-4
κₜ = DiscreteDiffusionFunction(κ, parameters=(κ₀ = 1e-2, ), loc=(Center, Center, Center))

PAR¹, PAR², PAR³ = Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid) #PAR bands
zₘₓₗ, zₑᵤ = Field{Center, Center, Center}(grid; indices = (:, :, 1:1)), Field{Center, Center, Center}(grid; indices = (:, :, 1:1)) #mixing and euphotic layer depth
Si̅ = Field{Center, Center, Center}(grid) #yearly maximum Si concentration
D_dust = Field{Center, Center, Center}(grid; indices = (:, :, 1:1)) #dust deposition rate

# Create a 'model' to run in Oceananignas
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                                                timestepper = :RungeKutta3,
                                                grid = grid,
                                                tracers = (:b, :T, :S, bgc.tracers...), #See https://clima.github.io/OceananigansDocumentation/stable/generated/ocean_wind_mixing_and_convection/ for T and S forcing
                                                coriolis = FPlane(f=1e-4),
                                                buoyancy = BuoyancyTracer(), 
                                                closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                                                forcing = bgc.forcing,
                                                boundary_conditions = bgc.boundary_conditions,
                                                auxiliary_fields = (; PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)
                                            )

@info "Model defined"

#=Pᵢ(x, y, z)= ((tanh((z+250)/100)+1)/2*(0.038)+0.002)*(12/14)*1e-6           #in molC L⁻¹
Zᵢ(x, y, z)= ((tanh((z+250)/100)+1)/2*(0.038)+0.008)*(12/14)*1e-6            #in molC L⁻¹
Chlᴾᵢ(x, y, z)= Pᵢ(x, y, z)*params.θᶜʰˡₘₐₓ.P/2       #in gChl L⁻¹
Chlᴰᵢ(x, y, z)= Pᵢ(x, y, z)*params.θᶜʰˡₘₐₓ.D/2       #in gChl L⁻¹
Feᴾᵢ(x, y, z)= Pᵢ(x, y, z)*params.θᶠᵉₒₚₜ.P       #in mol Fe L⁻¹
Feᴰᵢ(x, y, z)= Pᵢ(x, y, z)*params.θᶠᵉₒₚₜ.D       #in mol Fe L⁻¹
Siᴰᵢ(x, y, z)= Pᵢ(x, y, z)*params.θˢⁱᴰₘ      #in mol Si L⁻¹

#pretty arbitary
NO₃ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                  #in μmolN L⁻¹ (=mmol m⁻³)
NH₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*0.05+0.05             #in μmolN L⁻¹
PO₄ᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                  #in nmolP L⁻¹
Feᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                  #in nmolFe L⁻¹
Siᵢ(x, y, z)= (1-tanh((z+300)/150))/2*6+11.4                  #in mmolSi L⁻¹
Tᵢ(x, y, z) = 20 + 0.01*z
=#
set!(model, u = 0, v = 0, w = 0, b = 0, T = 17, S = 35, #Physical
                  P = 1e-9, D=1e-9, Z=1e-9, M=1e-9, Chlᴾ = 1e-9*params.θᶜʰˡₘₐₓ.P/2, Chlᴰ = 1e-9*params.θᶜʰˡₘₐₓ.D/2, Feᴾ = 1e-9*params.θᶠᵉₒₚₜ.P, Feᴰ = 1e-9*params.θᶠᵉₒₚₜ.D, Siᴰ = 1e-9*params.θˢⁱᴰₘ,
                  NO₃ = 5.0, NH₄ = 0.003, PO₄ = 5.0, Fe = 0.6, Si = 5.0, CaCO₃ = 0.0, 
                  DIC = 2200, Alk = 2400,#carbonates
                  O₂ = 240, Si̅ = 12, zₘₓₗ=100
    )

@info "Set initial values"

simulation = Simulation(model, Δt=1minute, stop_time=10years)

#=@inline H(t, t₀, t₁) = ifelse(t₀<t<t₁, 1.0, 0.0)
fmld1(t) = H.(t, 50days, 365days).*(1 ./(1 .+exp.(-(t-100days)/(5days)))).*(1 ./(1 .+exp.((t .-330days)./(25days))))

function analytical_zₘₓₗ(sim)
    #visual fit to north atlantic data
    t=sim.model.clock.time
    sim.model.auxiliary_fields.zₘₓₗ .= -(-10 -340 *(1 -fmld1(364.99999days)*exp(-t/25days)-fmld1(mod(t, 365days))))
end
simulation.callbacks[:update_zₘₓₗ] = Callback(analytical_zₘₓₗ, TimeInterval(0.5days))
=#
#PAR⁰(t) = 60*(1-cos((t+15days)*2π/(365days)))*(1 /(1 +0.2*exp(-((t-200days)/50days)^2))) .+ 2 #visual fit to north atlantic data, seems to be flattened in late summer/autumn so added exponential term, cloudier? Stdev much smaller when flatteing added (and both have mean error much smaller than stdev anyway)
PAR⁰(t) = 20
simulation.callbacks[:update_PAR] = Callback(Light.Morel.update, TimeInterval(1day), merge(Light.Morel.defaults, (; PAR⁰)))
simulation.callbacks[:nan_checker] = Callback(Oceananigans.Simulations.NaNChecker(merge(model.tracers, model.auxiliary_fields)))
Oceananigans.Simulations.erroring_NaNChecker!(simulation)
# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
    iteration(sim),
    prettytime(sim),
    prettytime(sim.Δt),
    prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10))

# Create a message to display every 100 timesteps                                
fields = Dict(zip(["$t" for t in [keys(model.tracers)..., keys(model.auxiliary_fields[(:PAR¹, :PAR², :PAR³)])...]], [model.tracers..., model.auxiliary_fields[(:PAR¹, :PAR², :PAR³)]...]))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, fields, filename="PISCES_example.nc", schedule=TimeInterval(.1days))
function check_neg(sim)
    for (name, tracer) in pairs((sim.model.tracers..., sim.model.auxiliary_fields...))
        if any(tracer .< 0)
            @warn "$name less than zero"
            #replace with kernal function 
            for k=1:sim.model.grid.Nz
                if tracer[1, 1, k]  < 0
                    tracer[:, : , k] .= 0.0
                end
            end
        end
    end  
end
simulation.callbacks[:negs] = Callback(check_neg)

simulation.output_writers[:diagnostics] = NetCDFOutputWriter(model, fields, filename="PISCES_example_diag.nc", schedule=SpecifiedTimes([0:1:100]))

@info "Simulation setup"

run!(simulation)
