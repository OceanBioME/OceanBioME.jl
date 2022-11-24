"
Integrates the biogeochemical models in a closed box
"
module BoxModel
using DiffEqBase, OrdinaryDiffEq, KernelAbstractions
using Oceananigans.Architectures: device

@kernel function calc_variable!(dy, y, params, t)
    p = @index(Global)
    yᵢ=[]
    for (i, tracer) in enumerate(params.variable_position)
        if tracer in params.forcing_dependencies[p]
            push!(yᵢ, y[i])
        end
    end
    dy[p] = params.forcing[p](0, 0, 0, t, yᵢ..., params.forcing_parameters.PAR(t), params.forcing_parameters)
end

function BoxModelRHS(dy, y, params, t)
    num = length(y)
    workgroup = min(num, 256)
    worksize = num

    calc_variable_kernal! = calc_variable!(device(params.architecture), workgroup, worksize)
    calc_variable_event = calc_variable_kernal!(dy, y, params, t)
    wait(calc_variable_event)
end

function run(model)
    return solve(model.problem, model.solver, dt=model.Δt, adaptive=model.adaptive)
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
=========
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