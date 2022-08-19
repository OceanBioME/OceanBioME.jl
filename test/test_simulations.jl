using OceanBioME, Test, Oceananigans, Printf
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

#forcing functions
t_function(x, y, z, t) = 15.0
s_function(x, y, z, t) = 35.0
surface_PAR(t) = 5

function run_simulation(bgc_model, optionsets, sediment, arch, c)
    grid = RectilinearGrid(arch, size=(1, 1, 10), extent=(1, 1, 100))
    params = getproperty(OceanBioME, bgc_model).defaults

    #setup optional tracers
    if optionsets == (:carbonates, )
        topboundaries = (DIC = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function)),)
    elseif optionsets == (:oxygen, )
        topboundaries = (OXY = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function)),)
    else
        topboundaries = NamedTuple()
    end

    #setup sediment
    if sediment
        sediment_bcs=Boundaries.setupsediment(grid)
    else
        sediment_bcs=NamedTuple()
    end

    #setup PAR 
    PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid) 

    #load bgc model
    bgc = Setup.Oceananigans(bgc_model, grid, params, PAR, optional_sets=optionsets, topboundaries=topboundaries, bottomboundaries=sediment_bcs, sinking=sediment) #with sinking and no sediment Oceananigans doesn't conserve the tracers (i.e. doesn't collect them at the bottom of the domain)

    #fix diffusivity
    κₜ(x, y, z, t) = 2e-2

    #setup model
    model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grid,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing =  bgc.forcing,
                            boundary_conditions = bgc.boundary_conditions,
                            auxiliary_fields =  !(sediment_bcs == NamedTuple()) ? merge((PAR=PAR, ), sediment_bcs.auxiliary_fields) : (PAR = PAR, )
    )

    #set default values for tracers
    for tracer in bgc.tracers
        set!(getproperty(model.tracers, tracer), getproperty(default_values, tracer))
    end
    
    set!(model, u=0, v=0, w=0, b=0)
    progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
    iteration(sim),
    prettytime(sim),
    prettytime(sim.Δt),
    prettytime(sim.run_wall_time))
    Δt = c*grid.Δzᵃᵃᶜ^2/κₜ(0, 0, 0, 0)
    duration = 100day

    simulation = Simulation(model, Δt=Δt, stop_time=duration) 
    simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))
    if sediment
        simulation.callbacks[:integrate_sediment] = sediment_bcs.callback
    end
    simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(10))

    simulation.output_writers[:profiles] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers, model.auxiliary_fields),
                          filename = "test/test_profiles.jld2",
                          indices = (1, 1, :),
                          schedule = TimeInterval(1days),     #TimeInterval(1days),
                          overwrite_existing = true)

    #get budgets and run simulation
    initial_budget = calculate_budget(model, sediment, getproperty(budget_tracers, bgc_model))
    run!(simulation)
    final_budget = calculate_budget(model, sediment, getproperty(budget_tracers, bgc_model))

    #test results loading
    results = OceanBioME.Plot.load_tracers(simulation)

    return bgc, model, simulation, initial_budget, final_budget, results
end

#for diagnostics
function record_budget!(sim, params)
    push!(budget, calculate_budget(sim.model, params.sediment, params.tracers))
end

function calculate_budget(model, sediment, tracers)
    budget = 0.0
    for tracer in tracers
        budget += sum(getproperty(model.tracers, tracer)[1, 1, 1:model.grid.Nz])
    end
    if sediment
        for sed_tracer in (:Nᵣ, :Nᵣᵣ, :Nᵣₑ)
            budget += getproperty(model.auxiliary_fields, sed_tracer)
        end
    end
    return budget
end

function calculate_budget(results::OceanBioME.Plot.model_results, bgc_model, sediment)
    budget = zeros(length(results.t))
    for tracer in getproperty(budget_tracers, bgc_model)
        ind = findfirst(results.tracers.=="$tracer")
        budget += sum(results.results[ind, 1, 1, :, :], dims=1)[1, :]
    end
    if sediment
        for sed_tracer in (:Nᵣ, :Nᵣᵣ, :Nᵣₑ)
            ind = findfirst(results.tracers.=="$tracer")
            budget += sum(results.results[ind, 1, 1, :, :], dims=1)[1, :]
        end
    end
    return budget
end

function calculate_C_budget(results::OceanBioME.Plot.model_results, bgc_model)
    budget = zeros(length(results.t))
    for tracer in getproperty(c_budget_tracers, bgc_model)
        ind = findfirst(results.tracers.=="$tracer")
        if !isnothing(ind) #incase carbonate model is turned off
            budget += sum(results.results[ind, 1, 1, 1:model.grid.Nz, :], dims=1)[1, :]
        end
    end
    if sediment
        sed_rd_sym = (:Rdᵣ, :Rdᵣᵣ, :Rd_red)
        for (j, sed_tracer) in enumerate((:Nᵣ, :Nᵣᵣ, :Nᵣₑ))
            ind = findfirst(results.tracers.=="$tracer")
            budget += sum(results.results[ind, 1, 1, 1:model.grid.Nz, :], dims=1)[1, :]*getproperty(OceanBioME.Boundaries.sediment, sed_rd_sym[j])
        end
    end

    airseaparams=merge(Boundaries.defaults.airseaflux, (gas=:CO₂, T=t_function, S=s_function))
    DIC_ind = findfirst(results.tracers.=="DIC")
    ALK_ind = findfirst(results.tracers.=="ALK")
    for (i, t) in enumerate(results.t)
         budget[i] += Boundaries.airseaflux(.5, .5, t, results.results[DIC_ind, 1, 1, end, i], results.results[ALK_ind, 1, 1, end, i], airseaparams) #if the example does not have a 1x1m x-y grid need to update this to be multipled by area
    end
    return budget
end

function run_test(bgc_model, optionsets, sediment, arch, c)
    bgc, model, simulation, ΣN₀, ΣN₁, results = run_simulation(bgc_model, optionsets, sediment, arch, c)
    @test ΣN₀≈ΣN₁
    for tracer in (model.tracers..., model.auxiliary_fields...)
        @test !any(isnan.(tracer[1, 1, 1:model.grid.Nz]))
        @test !any(isinf.(tracer[1, 1, 1:model.grid.Nz]))
        @test all(tracer[1, 1, 1:model.grid.Nz] .>= 0)
    end

    @test size(results.results)[1] == length(bgc.tracers) + length(model.auxiliary_fields)

    if :carbonates in optionsets
        c_budget = calculate_C_budget(results, bgc_model)
        @test c_budget[0] ≈ c_budget[1]
    end
end

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

models = (:LOBSTER, :NPZ)
arch = CPU()#GPU doesn't work for most at the moment but when it does itterate over that too
#global default values
Pᵢ(x, y, z)= 0.035
Zᵢ(x, y, z)= 0.045
Dᵢ(x, y, z)=0.0
DDᵢ(x, y, z)=0.0
NO₃ᵢ(x, y, z)= 12.0
NH₄ᵢ(x, y, z)= 0.05 
DOMᵢ(x, y, z)= 0 
DICᵢ(x, y, z)= 2200 
ALKᵢ(x, y, z)= 2400
OXYᵢ(x, y, z) = 240

const default_values = (P=Pᵢ, Z=Zᵢ, D=Dᵢ, DD=DDᵢ, NO₃=NO₃ᵢ, NH₄=NH₄ᵢ, DOM=DOMᵢ, DIC=DICᵢ, ALK=ALKᵢ, OXY=OXYᵢ)

budget_tracers = (LOBSTER = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), NPZ = (:N, :P, :Z))
c_budget_tracers = (LOBSTER = (P = LOBSTER.defaults.Rd_phy, Z=LOBSTER.defaults.Rd_phy, D=LOBSTER.defaults.Rd_phy, DD=LOBSTER.defaults.Rd_phy, DOM=LOBSTER.defaults.Rd_dom, DIC=1, ALK=1), )

c = 0.2
#=
@testset "Simulations" begin 
for model in models
    @info "Testing $model model"
    option_sets = keys(getproperty(OceanBioME, model).optional_tracers)
    for option in ([(o, ) for o in option_sets]..., ())
        for sediment in [true, false]
            if (sediment==false | (sediment==true && option == :oxygen))#sediment model depends on oxygen model
                @info "Testing optional tracers $option $(sediment ? "with" : "without") sediment"
                run_simulation(model, (option, ), sediment, arch, c)
            end
        end
    end
end
end=#