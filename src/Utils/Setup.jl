module Setup
export run!
using OceanBioME, Oceananigans, DiffEqBase, OrdinaryDiffEq, KernelAbstractions
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Architectures: device, arch_array

function loadmodel(model)
    models = propertynames(OceanBioME)
    deleteat!(models, findall(x->x in (:AirSeaFlux, :OceanBioME, :Light, :Particles), models))
    if !(model in models)
        throw(ArgumentError("Model '$model' is not an available model, choices are $models"))
    end

    return getproperty(OceanBioME, model)
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

    return (core=model.tracers, optional=optional_tracers)
end

@kernel function _setup_slip_vel!(vel, fv, z, params)
    k = @index(Global)
    @inbounds vel[:, :, k] .= fv(z[k], params)
end

function setuptracer(model, grid, tracer, field_dependencies, topboundaries, bottomboundaries, forcing_params; sinking, advection_scheme)
    forcing = Forcing(getproperty(model.forcing_functions, tracer), field_dependencies=field_dependencies, parameters=forcing_params)
    if (sinking && tracer in keys(model.sinking))
	slip_vel = arch_array(grid.architecture, zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2))
	workgroup = min(grid.Nz-1, 256)
        worksize = workgroup
        setup_slip_vel! = _setup_slip_vel!(device(grid.architecture), workgroup, worksize)
        setup_event = setup_slip_vel!(slip_vel, getproperty(model.sinking, tracer), grid.zᵃᵃᶜ, forcing_params)
        wait(setup_event)
        forcing = (forcing, AdvectiveForcing(advection_scheme(), w=slip_vel))
    end

    topboundary = tracer in keys(topboundaries) ? getproperty(topboundaries, tracer) : FluxBoundaryCondition(0)
    bottomboundary = tracer in keys(bottomboundaries) ? getproperty(bottomboundaries, tracer) : FluxBoundaryCondition(0)

    bcs = FieldBoundaryConditions(top = topboundary, bottom = bottomboundary)
    
    return forcing, bcs
end

function Oceananigans(model::Symbol, 
                                    grid::AbstractGrid,  
                                    params::NamedTuple, 
                                    PAR; #should add a test for this to be an auxiliary field or a fucntion (but can't distinguish an aux field or tracer field)
                                    topboundaries::NamedTuple=NamedTuple(),
                                    bottomboundaries::NamedTuple=NamedTuple(),
                                    optional_sets::Tuple=(),
                                    sinking = true, 
                                    advection_scheme = UpwindBiasedFifthOrder)

    model=loadmodel(model)
    tracers=loadtracers(model, optional_sets)

    forcing_functions=()
    boundary_functions=()

    forcing_params=merge(params, (PAR = PAR, ))

    for tracer in tracers.core
        forcing, bcs = setuptracer(model, grid, tracer, tracers.core, topboundaries, bottomboundaries, forcing_params; sinking=sinking, advection_scheme=advection_scheme)
        forcing_functions = (forcing_functions..., forcing)
        boundary_functions = (boundary_functions..., bcs)
    end

    for optionset in tracers.optional
        for tracer in optionset
            forcing, bcs = setuptracer(model, grid, tracer, (tracers.core..., optionset...), topboundaries, bottomboundaries, forcing_params; sinking=sinking, advection_scheme=advection_scheme)
            forcing_functions = (forcing_functions..., forcing)
            boundary_functions = (boundary_functions..., bcs)
        end
    end

    tracers = (tracers.core..., (tracers.optional...)...)
    forcing = (; zip(tracers, forcing_functions)...)
    boundaries = (; zip(tracers, boundary_functions)...)

    return (tracers=tracers, forcing=forcing, boundary_conditions=boundaries)
end

function BoxModel(model::Symbol, params::NamedTuple, initial_values::NamedTuple, tᵢ::AbstractFloat, tₑ::AbstractFloat; optional_sets::Tuple=(), Δt=200, solver=Euler(), adaptive=false, architecture=CPU())
    model=loadmodel(model)
    tracers=loadtracers(model, optional_sets)

    if length(initial_values) != length((tracers.core..., (tracers.optional...)...))
        throw(ArgumentError("The incorrect number of initial values have been provided, need initial values for $tracers"))
    end

    y₀ = []
    forcing_functions=[]
    forcing_dependencies=[]
    variable_position=[]
    for tracer in tracers.core
        push!(forcing_functions, getproperty(model.forcing_functions, tracer))
        push!(y₀, getproperty(initial_values, tracer))
        push!(forcing_dependencies, tracers.core)
        push!(variable_position, tracer)
    end

    for optionset in tracers.optional
        for tracer in optionset
            push!(forcing_functions, getproperty(model.forcing_functions, tracer))
            push!(y₀, getproperty(initial_values, tracer))
            push!(forcing_dependencies, (tracers.core..., optionset...))
            push!(variable_position, tracer)
        end
    end
    params = (forcing=forcing_functions, forcing_dependencies = forcing_dependencies, forcing_parameters=params, variable_position=variable_position, architecture=architecture)
    problem = ODEProblem(OceanBioME.BoxModel.BoxModelRHS, vcat(y₀...), (tᵢ, tₑ), params)
    return (tracers=(tracers.core..., (tracers.optional...)...), parameters = params, forcing=forcing_functions, problem=problem, Δt=Δt, solver=solver, adaptive=adaptive, solution=vcat([]...), t=vcat([]...))
end
end
