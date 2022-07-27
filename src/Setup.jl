module Setup
using BGC, Oceananigans
using Oceananigans.Grids: AbstractGrid

function Oceananigans(model::Symbol, 
                                    grid::AbstractGrid,  
                                    params::NamedTuple, 
                                    PAR; #should add a test for this to be an auxiliary field or a fucntion (but can't distinguish an aux field or tracer field)
                                    topboundaries::NamedTuple=NamedTuple(),
                                    bottomboundaries::NamedTuple=NamedTuple(),
                                    optional_sets::Tuple=())

    models = propertynames(BGC)
    deleteat!(models, findall(x->x in (:AirSeaFlux, :BGC, :Light, :Particles), models))
    if !(model in models)
        throw(ArgumentError("Model '$model' is not an available model, choices are $models"))
    end

    model = getproperty(BGC, model)

    optional_tracers = []
    available_optionsets = keys(model.optional_tracers)
    for option in optional_sets
        if !(option in available_optionsets)
            throw(ArgumentError("$option is not an available option set in the model '$model', available are $available_optionsets"))
        end
        push!(optional_tracers, getproperty(model.optional_tracers, option))
    end

    tracers=(model.tracers..., (optional_tracers...)...)

    forcing_functions=()
    boundary_functions=()

    forcing_params=merge(params, (PAR = PAR, ))

    for tracer in tracers
        forcing = Forcing(getproperty(model.forcing_functions, tracer), field_dependencies=tracers, parameters=forcing_params)
        if tracer in keys(model.sinking)
            forcing = (forcing, AdvectiveForcing(WENO5(; grid), w=getproperty(model.sinking, tracer)))
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
end