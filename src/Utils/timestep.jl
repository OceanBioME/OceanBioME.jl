using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Grids: Center, znodes, zspacing, minimum_zspacing

@inline column_advection_timescale(model) = minimum_zspacing(model.grid) / maximum_sinking_velocity(model.biogeochemistry)

@inline sinking_advection_timescale(model) = min(minimum_zspacing(model.grid) / maximum_sinking_velocity(model.biogeochemistry), cell_advection_timescale(model))
