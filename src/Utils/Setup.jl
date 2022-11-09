module Setup
export run!
using OceanBioME, Oceananigans, DiffEqBase, OrdinaryDiffEq
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Architectures: arch_array
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition, ValueBoundaryCondition
using Oceananigans.Fields: ZFaceField

#####
##### Infastructure to load biogeochemical models
#####

function loadmodel(model)
    models = propertynames(OceanBioME)
    deleteat!(models, findall(x->x in (:AirSeaFlux, :OceanBioME, :Light, :Particles), models))
    if !(model in models)
        throw(ArgumentError("Model '$model' is not an available model, choices are $models"))
    end
    properties = getproperty(OceanBioME, model)

    #horrible but don't know how else todo this
    discrete_forcings = NamedTuple()
    try 
        discrete_forcings = properties.discrete_forcing
    catch
        true
    end

    sinking = NamedTuple()
    try 
        sinking = properties.sinking
    catch
        true
    end

    optional = NamedTuple()
    try 
        optional = properties.optional_tracers
    catch
        true
    end

    dependencies = properties.tracers
    try
        dependencies = (dependencies..., properties.required_fields...)
    catch
        true
    end
    return properties, dependencies, discrete_forcings, dependencies, sinking, optional
end

function loadtracers(model, optional_tracers, optional_sets)
    optionals = []
    for option in optional_sets
        if !(option in keys(optional_tracers))
            throw(ArgumentError("$option is not an available option set, available are $(keys(optional_tracers)))"))
        end
        push!(optionals, getproperty(optional_tracers, option))
    end

    return (core=model.tracers, optional=optionals)
end

@inline get_sinking_vel(slip_vel::Function, z, args...) = slip_vel(z, args...)
@inline get_sinking_vel(slip_vel::Number, z, args...) = slip_vel

function setuptracer(model, grid, tracer, field_dependencies, topboundaries, bottomboundaries, forcing_params, discrete_forcings, sinking_forcings; sinking, open_bottom)
    discrete = ifelse(tracer in keys(discrete_forcings), getproperty(discrete_forcings, tracer), false)
    forcing = Forcing(getproperty(model.forcing_functions, tracer), field_dependencies=field_dependencies, parameters=forcing_params, discrete_form = discrete)
    if (sinking && tracer in keys(sinking_forcings))
        w_slip = ZFaceField(grid)
        for k=1:grid.Nz 
            @inbounds w_slip[:, :, k] .= get_sinking_vel(model.sinking[tracer], grid.zᵃᵃᶠ[k], forcing_params)*(open_bottom ? 1.0 : (1 - exp((1-k)/2)))
        end
        forcing = (forcing, AdvectiveForcing(WENO(;grid), w=w_slip))
    end

    topboundary = tracer in keys(topboundaries) ? getproperty(topboundaries, tracer) : FluxBoundaryCondition(0)
    bottomboundary = tracer in keys(bottomboundaries) ? getproperty(bottomboundaries, tracer) : FluxBoundaryCondition(0)

    bcs = FieldBoundaryConditions(top = topboundary, bottom = bottomboundary)
    
    return forcing, bcs
end

"""
    Setup.Oceananigans(model::Symbol, 
                                    grid::AbstractGrid,  
                                    forcing_params::NamedTuple; 
                                    topboundaries::NamedTuple=NamedTuple(),
                                    bottomboundaries::NamedTuple=NamedTuple(),
                                    optional_sets::Tuple=(),
                                    sinking = true, 
                                    supress_required_fields_warning = false,
                                    open_bottom=true)

Returns the tracers, forcing functions, and boundary conditions for a biogeochemical model.
Arguments
==========
* `model`: Symbol name of model to be used (e.g. `:LOBSTER`)
* `grid`: Oceananigans grid (used to setup sinking if present)
* `forcing_params`: Parameters for model, usually `model.defaults` (no way for us to automagically get that for now)
Keyword arguments
==============
* `topboundaries` and `bottomboundaries`: NamedTuple of boundary conditions (values) for tracers (keys)
* `optional_sets`: Tuple of Symbols with names of optional tracer sets, you can usually see what is available by `MODEL_NAME.optional_tracers`
* `sinking`: turn off particle sinking
* `supress_required_fields_warning`: turn off the warning of additional required model fields (e.g. PAR)
* `open_bottom`: if true particles may fall out of model, if false then sinking speed is smoothly brought to zero near bottom boundary (may cause odd behaviour/very slow sinking on small grids)
"""
function Oceananigans(model::Symbol, 
                                    grid::AbstractGrid,  
                                    forcing_params::NamedTuple; 
                                    topboundaries::NamedTuple=NamedTuple(),
                                    bottomboundaries::NamedTuple=NamedTuple(),
                                    optional_sets::Tuple=(),
                                    sinking = true, 
                                    supress_required_fields_warning = false,
                                    open_bottom=true)

    model, dependencies, discrete_forcings, dependencies, sinking_forcings, optional_tracers = loadmodel(model)
    tracers=loadtracers(model, optional_tracers, optional_sets)

    discrete_forcings = ifelse(length(discrete_forcings) == 0, NamedTuple{(tracers.core..., (tracers.optional...)...)}(repeat([false], length((tracers.core..., (tracers.optional...)...)))), discrete_forcings)

    forcing_functions = ()
    boundary_functions = ()

    for tracer in tracers.core
        forcing, bcs = setuptracer(model, grid, tracer, dependencies, topboundaries, bottomboundaries, forcing_params, discrete_forcings, sinking_forcings; sinking=sinking, open_bottom=open_bottom)
        forcing_functions = (forcing_functions..., forcing)
        boundary_functions = (boundary_functions..., bcs)
    end

    for optionset in tracers.optional
        for tracer in optionset
            forcing, bcs = setuptracer(model, grid, tracer, (dependencies..., optionset...), topboundaries, bottomboundaries, forcing_params, discrete_forcings, sinking_forcings; sinking=sinking, open_bottom=open_bottom)
            forcing_functions = (forcing_functions..., forcing)
            boundary_functions = (boundary_functions..., bcs)
        end
    end

    tracers = (tracers.core..., (tracers.optional...)...)
    forcing = (; zip(tracers, forcing_functions)...)
    boundaries = (; zip(tracers, boundary_functions)...)
    
    if !supress_required_fields_warning
        try
            @warn "This model requires $(model.required_fields) to be separatly defined (as tracer or auxiliary fields)"
        catch
            true
        end

        try
            @warn "This model requires $(model.requried_parameters) to be separatly defined in addition to the default parameters (MODEL_NAME.defaults)"
        catch
            true
        end
    end
    return (tracers=tracers, forcing=forcing, boundary_conditions=boundaries)
end

"""
    Setup.Oceananigans(model::Symbol, 
                                    params::NamedTuple, 
                                    initial_values::NamedTuple, 
                                    tᵢ::AbstractFloat, 
                                    tₑ::AbstractFloat; 
                                    optional_sets::Tuple=(), 
                                    Δt=200, 
                                    solver=Euler(), 
                                    adaptive=false, 
                                    architecture=CPU())

Returns a box model (actual type coming soon) which can be solved by `BoxModel.run`.
Arguments
==========
* `model`: Symbol name of model to be used (e.g. `:LOBSTER`)
* `params`: Parameters for model, usually `model.defaults` (no way for us to automagically get that for now)
* `initial_vales`: NamedTuple of initial values for tracers
* `tᵢ`: start time
* `tₑ`: end time
Keyword arguments
==============
* `optional_sets`: Tuple of Symbols with names of optional tracer sets, you can usually see what is available by `MODEL_NAME.optional_tracers`
* `Δt`: time step for integrator
* `solver`: `OrdinaryDiffEq` solver 
* `adaptive`: adapt timestep length?
* `architecture`: `KernelAbstractions` architecture for distribution of computation
"""

function BoxModel(model::Symbol, 
                              params::NamedTuple, 
                              initial_values::NamedTuple, 
                              tᵢ::AbstractFloat, 
                              tₑ::AbstractFloat; 
                              optional_sets::Tuple=(), 
                              Δt=200, 
                              solver=Euler(), 
                              adaptive=false, 
                              architecture=CPU())

    model, dependencies, discrete_forcings, sinking_forcings, optional_tracers=loadmodel(model)
    tracers=loadtracers(model, optional_tracers, optional_sets)

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
