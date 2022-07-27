module Setup
export run!
using BGC, Oceananigans, DiffEqBase, OrdinaryDiffEq
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Architectures: arch_array

function loadmodel(model)
    models = propertynames(BGC)
    deleteat!(models, findall(x->x in (:AirSeaFlux, :BGC, :Light, :Particles), models))
    if !(model in models)
        throw(ArgumentError("Model '$model' is not an available model, choices are $models"))
    end

    return getproperty(BGC, model)
end

function loadtracers(model, optional_sets)
    optional_tracers = []
    available_optionsets = keys(model.optional_tracers)
    for option in optional_sets
        if !(option in available_optionsets)
            throw(ArgumentError("$option is not an available option set in the model '$model', available are $available_optionsets"))
        end
        push!(optional_tracers, getproperty(model.optional_tracers, option))
    end

    return (model.tracers..., (optional_tracers...)...)
end

function Oceananigans(model::Symbol, 
                                    grid::AbstractGrid,  
                                    params::NamedTuple, 
                                    PAR; #should add a test for this to be an auxiliary field or a fucntion (but can't distinguish an aux field or tracer field)
                                    topboundaries::NamedTuple=NamedTuple(),
                                    bottomboundaries::NamedTuple=NamedTuple(),
                                    optional_sets::Tuple=())

    model=loadmodel(model)
    tracers=loadtracers(model, optional_sets)

    forcing_functions=()
    boundary_functions=()

    forcing_params=merge(params, (PAR = PAR, ))

    for tracer in tracers
        forcing = Forcing(getproperty(model.forcing_functions, tracer), field_dependencies=tracers, parameters=forcing_params)
        if tracer in keys(model.sinking)
            slip_vel = arch_array(grid.architecture, zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2))
            for k=0:grid.Nz-2
                @inbounds slip_vel[:, :, k] .= getproperty(model.sinking, tracer)(grid.zᵃᵃᶜ[k], params)
            end
            forcing = (forcing, AdvectiveForcing(WENO5(; grid), w=slip_vel))
        end
        forcing_functions = (forcing_functions..., forcing)

        topboundary = tracer in keys(topboundaries) ? getproperty(topboundaries, tracer) : FluxBoundaryCondition(0)
        bottomboundary = tracer in keys(bottomboundaries) ? getproperty(bottomboundaries, tracer) : FluxBoundaryCondition(0)

        bcs = FieldBoundaryConditions(top = topboundary, bottom = bottomboundary)

        boundary_functions = (boundary_functions..., bcs)
    end

    forcing = (; zip(tracers, forcing_functions)...)
    boundaries = (; zip(tracers, boundary_functions)...)

    return (tracers=tracers, forcing=forcing, boundary_conditions=boundaries)
end

mutable struct BoxModel
    tracers
    forcing
    problem
    Δt
    solver
    adaptive
    solution
end

function BoxModelRHS(y, params, t)
    return (vcat([params.forcing[i](0, 0, 0, t, y..., params) for i in 1:length(y)]...))
end

function BoxModel(model::Symbol, params::NamedTuple, initial_values::NamedTuple, tᵢ::AbstractFloat, tₑ::AbstractFloat; optional_sets::Tuple=(), Δt=200, solver=RK4(), adaptive=false)
    model=loadmodel(model)
    tracers=loadtracers(model, optional_sets)

    if length(initial_values) != length(tracers)
        throw(ArgumentError("The incorrect number of initial values have been provided, need initial values for $tracers"))
    end
    y₀ = values(initial_values)

    forcing_functions=()
    for tracer in tracers
        forcing = getproperty(model.forcing_functions, tracer)
    end

    params = merge(params, (forcing=forcing_functions, ))
    problem = ODEProblem(BoxModelRHS, y₀, (tᵢ, tₑ), params)
    return BoxModel(tracers, forcing, problem, Δt, solver, adaptive, [])
end

function run!(model::BoxModel)
    solution = solve(model.problem, model.solver, dt=model.Δt, adaptive=model.adaptive)
    model.solution = solution
end
end