module CarbonateSystemModels

export SimpleCarbonateSystem

@inline carbonte_tendencies(cs, bgc, i, j, k, grid, val_name, clock, fields, auxiliary_fields) = zero(grid)

include("simple_carbonate.jl")

end