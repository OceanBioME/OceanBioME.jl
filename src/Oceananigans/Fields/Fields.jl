module Fields

using Oceananigans.Fields: ZeroField, Field

import Statistics: mean

include("field_tuples.jl")
include("interpolate.jl")

mean(f::ZeroField; kwargs...) = 0.0

@inline isacolumn(field::Field) =  size(field)[3] == 1

end # module