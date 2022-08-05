"""
LOBSTER model based on
2005 A four-dimensional mesoscale map of the spring bloom in the northeast Atlantic (POMME experiment): Results of a prognostic model
2001 Impact of sub-mesoscale physics on production and subduction of phytoplankton in an oligotrophic regime
2012How does dynamical spatial variability impact 234Th-derived estimates of organic export
using flux boundary condition
ignore aggregation term 
annual cycle 
add callback to diagnose pCO2
"""

using Random
using Printf
using Plots
using JLD2
using NetCDF
using HDF5
using Interpolations
using Statistics

using Oceananigans#9e8cae18-63c1-5223-a75c-80ca9d6e9a09
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

using BGC, SugarKelp, StructArrays

params = LOBSTER.defaults


begin#Load data
    ######## 1: import annual cycle physics data  #Global Ocean 1/12° Physics Analysis and Forecast updated Daily
    filename1 = "subpolar_physics.nc" # subtropical_physics.nc  subpolar_physics.nc
    #ncinfo(filename1)
    time_series_second = (0:364)days # start from zero if we don't use extrapolation, we cannot use extrapolation if we wana use annual cycle  
    so = ncread(filename1, "so");  #salinity
    so_scale_factor = ncgetatt(filename1, "so", "scale_factor") #Real_Value = (Display_Value X scale_factor) + add_offset
    so_add_offset = ncgetatt(filename1, "so", "add_offset")
    salinity = mean(so, dims=(1,2))[1:365]*so_scale_factor.+so_add_offset # use [1:365] cause sometimes one year has 366 days. 
    #converted to interpolations to access them at arbitary time, how to use it: salinity_itp(mod(timeinseconds,364days))
    #plot(salinity)
    thetao = ncread(filename1, "thetao");  #temperature
    thetao_scale_factor = ncgetatt(filename1, "thetao", "scale_factor") 
    thetao_add_offset = ncgetatt(filename1, "thetao", "add_offset")
    temperature = mean(thetao, dims=(1,2))[1:365]*thetao_scale_factor.+thetao_add_offset
    #plot(temperature)
    mlotst = ncread(filename1, "mlotst"); #mixed_layer_depth
    mlotst_scale_factor = ncgetatt(filename1, "mlotst", "scale_factor") 
    mlotst_add_offset = ncgetatt(filename1, "mlotst", "add_offset")
    mixed_layer_depth = mean(mlotst, dims=(1,2))[1:365]*mlotst_scale_factor.+mlotst_add_offset
    #plot(mixed_layer_depth)

    ######## 2: import annual cycle chl data  #Global Ocean Biogeochemistry Analysis and Forecast
    filename2 = "subpolar_chl.nc"    #subpolar_chl.nc
    chl = ncread(filename2, "chl");  #chl scale_factor=1 add_offset=089639299014
    chl_mean = mean(chl, dims=(1,2))[1,1,:,1:365] # mg m-3, unit no need to change. 
    depth_chl = ncread(filename2, "depth");
    #heatmap(1:365, -depth_chl[end:-1:1], chl_mean[end:-1:1,:])

    ######## 3: import annual cycle PAR data #Ocean Color  VIIRS-SNPP PAR daily 9km
    path="./subpolar/"    #subtropical   #./subpolar/
    par_mean_timeseries=zeros(365)
    for i in 1:365    #https://discourse.julialang.org/t/leading-zeros/30450
        string_i = lpad(string(i), 3, '0')
        filename3=path*"V2020"*string_i*".L3b_DAY_SNPP_PAR.x.nc"
        fid = h5open(filename3, "r")
        par=read(fid["level-3_binned_data/par"])
        BinList=read(fid["level-3_binned_data/BinList"])  #(:bin_num, :nobs, :nscenes, :weights, :time_rec) 
        par_mean_timeseries[i] = mean([par[i][1]/BinList[i][4] for i in 1:length(par)])*3.99e-10*545e12/(1day)  #from einstin/m^2/day to W/m^2
    end
end
salinity_itp = LinearInterpolation(time_series_second, salinity) 
temperature_itp = LinearInterpolation(time_series_second, temperature)  
surface_PAR_itp = LinearInterpolation((0:364)day, par_mean_timeseries)
mld_itp = LinearInterpolation(time_series_second, mixed_layer_depth)  
surface_PAR(t) = surface_PAR_itp(mod(t, 364days))

t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))

# Simulation duration  
duration=2years
# Define the grid

Lx = 1   #500
Ly = 500
Nx = 1
Ny = 1
Nz = 150
Lz = 600 # domain depth             # subpolar mixed layer depth max 427m 

#inspired by the grid spacing in Copurnicus' models
#spaced_z_faces(k) = 0.5*(1-exp((7.1388/Nz)*(Nz-k)))
#=
#inspired by example
stretching = 5
refinement = 2.5
h(k) = (k-1)/ Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(
                size=(Nx, Ny, Nz), 
                x=(0,Lx),
                y=(0,Ly),
                z=z_faces)=#

grid = RectilinearGrid(size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz))

PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)
begin #setup bouyancy
    B₀ = 0e-8    #m²s⁻³  B₀ = 4.24e-8  
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
end

#Load the BGC model
dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
bgc = Setup.Oceananigans(:LOBSTER, grid, params, PAR, topboundaries=(DIC=dic_bc, ), optional_sets=(:carbonates, ))

@info "Setup BGC model"

begin #Setup the kelp particles
    ##refactor equations so they return in the correct units etc
    function sugarkelpequations(x, y, z, t, A, N, C, NO₃, irr, params, Δt)
        if !iszero(A)
            temp=params.temp(mod(t, 364days))

            p=SugarKelp.p(temp, irr, params.paramset)
            e=SugarKelp.e(C, params.paramset)
            μ = SugarKelp.eval_μ(A, N, C, temp, params.λ[1+floor(Int, mod(t, 364days)/day)], params.paramset)
            j = SugarKelp.eval_j(NO₃, N, params.urel, params.paramset)
            ν = SugarKelp.ν(A, params.paramset)
            dA = (μ - ν) * A / (60*60*24)
            dN = (j / params.paramset.K_A - μ * (N + params.paramset.N_struct)) / (60*60*24)
            dC = ((p* (1 - e) - SugarKelp.r(temp, μ, j, params.resp_model, params.paramset)) / params.paramset.K_A - μ * (C + params.paramset.C_struct)) / (60*60*24)

            A_new = A+dA*Δt 
            N_new = N+dN*Δt 
            C_new = C+dC*Δt

            if C_new < params.paramset.C_min
                A_new *= (1-(params.paramset.C_min - C)/params.paramset.C_struct)
                C_new = params.paramset.C_min
            end

            p *= A_new / (60*60*24*12*0.001)#gC/dm^2/hr to mmol C/s
            e *= p#mmol C/s
            ν *= A_new*(N_new + params.paramset.N_struct) / (60*60*24*14*0.001)#1/hr to mmol N/s
            j *= A_new / (60*60*24*14*0.001)#gN/dm^2/hr to mmol N/s
        else
            A_new, N_new, C_new, j, p, e, ν = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        end
        return (A = A_new, N = N_new, C = C_new, j = j, p = p, e = e, ν = ν)
    end

    lat=57.5
    λ_arr=SugarKelp.gen_λ(lat)

    struct Kelp
        #position
        x :: AbstractFloat
        y :: AbstractFloat
        z :: AbstractFloat
        #velocity
        u :: AbstractFloat
        v :: AbstractFloat
        w :: AbstractFloat
        #properties
        ##integral fields
        A :: AbstractFloat
        N :: AbstractFloat
        C :: AbstractFloat
        ##diagnostic fields
        j :: AbstractFloat
        p :: AbstractFloat
        e :: AbstractFloat
        ν :: AbstractFloat
        #tracked fields
        NO₃  :: AbstractFloat
        PAR :: AbstractFloat
    end

    n_kelp = 100

    x₀ₖ = (Lx/2)*ones(n_kelp)
    y₀ₖ = (Ly/2)*ones(n_kelp)
    z₀ₖ = [-100:-1;]

    #incluidng velocity for later use
    u₀ₖ = zeros(n_kelp)
    v₀ₖ = zeros(n_kelp)
    w₀ₖ = zeros(n_kelp)

    #start it off dead so it doesn't effect the BGC burnin
    a₀ = zeros(n_kelp)
    n₀ = zeros(n_kelp)
    c₀ = zeros(n_kelp)
    j₀ = zeros(n_kelp)
    p₀ = zeros(n_kelp)
    e₀ = zeros(n_kelp)
    ν₀ = zeros(n_kelp)
    NO₃₀ = zeros(n_kelp)
    PAR₀ = zeros(n_kelp)

    particles = StructArray{Kelp}((x₀ₖ, y₀ₖ, z₀ₖ, u₀ₖ, v₀ₖ, w₀ₖ, a₀, n₀, c₀, j₀, p₀, e₀, ν₀, NO₃₀, PAR₀))

    source_fields = ((tracer=:NO₃, property=:NO₃, scalefactor=1.0), 
                            (tracer=:PAR, property=:PAR, scalefactor=1.0))
    sink_fields = ((tracer=:NO₃, property=:j, scalefactor=-1.0), 
                        (tracer=:DIC, property=:p, scalefactor=-1.0),
                        (tracer=:DOM, property=:e, scalefactor=1.0/params.Rd_dom),
                        (tracer=:DD, property=:ν, scalefactor=1.0))
    #need to change urel at some point
    kelp_particles = Particles.setup(particles, sugarkelpequations, 
                                        (:A, :N, :C, :NO₃, :PAR), #forcing function property dependencies
                                        (urel=0.15, temp=temperature_itp, λ=λ_arr, resp_model=2, paramset=SugarKelp.broch2013params), #forcing function parameters
                                        (), #forcing function integrals - changed all kelp model variables to diagnostic since it is much easier to enforce the extreme carbon limit correctly
                                        (:A, :N, :C, :j, :p, :e, :ν), #forcing function diagnostic fields
                                        source_fields,
                                        sink_fields,
                                        100.0)#density (100 per meter down)
end
@info "Defined kelp particles"
#κₜ(x, y, z, t) = 1e-2*max(1-(z+50)^2/50^2,0)+1e-5;
@inline κₜ(x, y, z, t) = 1e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-5;

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
                            auxiliary_fields = (PAR=PAR, ),
                            particles = kelp_particles)#comment out this line if using functional form of PAR
@info "Setup model"
#set initial conditions
initial_mixed_layer_depth = -100 # m
stratification(z) = z < initial_mixed_layer_depth ? N² * z : N² * (initial_mixed_layer_depth)
bᵢ(x, y, z) = stratification(z)         #+ 1e-1 * Ξ(z) * N² * model.grid.Lz

#could initiate P from chl data?
Pᵢ(x,y,z)= (tanh((z+250)/100)+1)/2*(0.038)+0.002          # ((tanh((z+100)/50)-1)/2*0.23+0.23)*16/106  
Zᵢ(x,y,z)= (tanh((z+250)/100)+1)/2*(0.038)+0.008          # ((tanh((z+100)/50)-1)/2*0.3+0.3)*16/106         
Dᵢ(x,y,z)=0
DDᵢ(x,y,z)=0
NO₃ᵢ(x,y,z)= (1-tanh((z+300)/150))/2*6+11.4   #  # 17.5*(1-tanh((z+100)/10))/2
NH₄ᵢ(x,y,z)= (1-tanh((z+300)/150))/2*0.05+0.05       #1e-1*(1-tanh((z+100)/10))/2
DOMᵢ(x,y,z)= 0 
DICᵢ(x,y,z)= 2380   #  mmol/m^-3
ALKᵢ(x,y,z)= 2720   #  mmol/m^-3

## `set!` the `model` fields using functions or constants:
set!(model, b=bᵢ, P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ,DIC=DICᵢ,ALK=ALKᵢ, u=0, v=0, w=0)

# ## Setting up a simulation
simulation = Simulation(model, Δt=200, stop_time=6*31day+3*30day+28*day)  #Δt=0.5*(Lz/Nz)^2/1e-2,
#run to start of november

simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)));#comment out if using PAR functiuon

## Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

# Vertical slice
simulation.output_writers[:profiles] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, model.auxiliary_fields),
                          filename = "profile_subpolar.jld2",
                          indices = (1, 1, :),
                          schedule = TimeInterval(1days),     #TimeInterval(1days),
                            overwrite_existing = true)

#checkpoint after burnin so we can just reset to there next time
simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([6*31day+3*30day+28*day]), prefix="model_checkpoint")
@info "Setup simulation"

run!(simulation)
#reset particles

model.particles.properties.A .= 30.0*ones(n_kelp)
model.particles.properties.N .= 0.01*ones(n_kelp)
model.particles.properties.C .= 0.1*ones(n_kelp)

#start recording particles
simulation.output_writers[:particles] = JLD2OutputWriter(model, (particles=model.particles,), 
                          filename = "particles.jld2",
                          schedule = TimeInterval(1day),
                          overwrite_existing = true)

simulation.stop_time = duration

run!(simulation)
@info "Simulation finished"
## Coordinate arrays
xw, yw, zw = nodes(model.velocities.w)
xb, yb, zb = nodes(model.tracers.b)

results = BGC.Plot.load_tracers(simulation)
BGC.Plot.profiles(results)
savefig("annual_cycle_subpolar_highinit_kelp.pdf")

particles = BGC.Plot.load_particles(simulation)
BGC.Plot.particles(particles)
savefig("annual_cycle_subpolar_highinit_kelp_particles.pdf")