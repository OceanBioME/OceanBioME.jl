using Oceananigans.Fields: TendencyFields

import Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper

function QuasiAdamsBashforth2TimeStepper(grid, tracers, forced_auxiliary_fields,
        χ = 0.1;
        implicit_solver::IT = nothing,
        Gⁿ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid)),
        G⁻ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid))) where IT

    FT = eltype(grid)
    GT = typeof(Gⁿ)

    return QuasiAdamsBashforth2TimeStepper{FT, GT, IT}(χ, Inf, Gⁿ, G⁻, implicit_solver)
end