module Models

export Sediments, NPZD, NutrientPhytoplanktonZooplanktonDetritus, LOBSTER, SLatissima, GasExchange

include("Sediments/Sediments.jl")
include("AdvectedPopulations/LOBSTER/LOBSTER.jl")
include("AdvectedPopulations/NPZD.jl")
include("Individuals/SLatissima.jl")
include("gasexchange.jl")

using .Sediments
using .LOBSTERModel
using .NPZDModel
using .SLatissimaModel

end # module