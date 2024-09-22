include("dependencies_for_runtests.jl")

include("test_utils.jl")
include("test_light.jl")
include("test_slatissima.jl")
include("test_LOBSTER.jl")
include("test_NPZD.jl")
include("test_PISCES.jl")
include("test_gasexchange_carbon_chem.jl")
include("test_sediments.jl")

if architecture == CPU() 
    include("test_boxmodel.jl") # box models (probably) don't work on GPU, and it wouldn't be faster anyway
end

architecture == CPU() && @testset "Doctests" begin
    doctest(OceanBioME)
end
