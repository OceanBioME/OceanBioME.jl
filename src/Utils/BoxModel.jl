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

end