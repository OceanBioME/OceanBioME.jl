using Oceananigans.Fields: TendencyFields

import Oceananigans.TimeSteppers: RungeKutta3TimeStepper

function RungeKutta3TimeStepper(grid, tracers, forced_auxiliary_fields;
    implicit_solver::TI = nothing,
    Gⁿ::TG = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid)),
    G⁻ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid))) where {TI, TG}

    !isnothing(implicit_solver) &&
    @warn("Implicit-explicit time-stepping with RungeKutta3TimeStepper is not tested. " * 
    "\n implicit_solver: $(typeof(implicit_solver))")

    γ¹ = 8 // 15
    γ² = 5 // 12
    γ³ = 3 // 4

    ζ² = -17 // 60
    ζ³ = -5 // 12

    FT = eltype(grid)

    return RungeKutta3TimeStepper{FT, TG, TI}(γ¹, γ², γ³, ζ², ζ³, Gⁿ, G⁻, implicit_solver)
end