using OceanBioME, Documenter, Test, CUDA, Oceananigans

architecture = CUDA.has_cuda() ? GPU() : CPU()

include("test_utils.jl")
include("test_light.jl")
include("test_slatissima.jl")
include("test_LOBSTER.jl")
include("test_NPZD.jl")
include("test_gasexchange.jl")
include("test_sediments.jl")

architecture == CPU() && @testset "Doctests" begin
    doctest(OceanBioME)
end
