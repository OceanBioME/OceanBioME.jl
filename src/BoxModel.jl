module BoxModel
using DiffEqBase, OrdinaryDiffEq
function BoxModelRHS(dy, y, params, t)
    dy .= (vcat([params.forcing[i](0, 0, 0, t, y..., params.parameters) for i in 1:length(y)]...))
end

function run(model)
    return solve(model.problem, model.solver, dt=model.Δt, adaptive=model.adaptive)
end
end