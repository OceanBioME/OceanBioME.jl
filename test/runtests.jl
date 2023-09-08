using OceanBioME, Documenter, Test, Oceananigans

if !(@isdefined arch)
    arch = CPU()
end

include("test_utils.jl")
include("test_light.jl")
include("test_LOBSTER.jl")
include("test_NPZD.jl")
include("test_gasexchange.jl")
include("test_slatissima.jl")
include("test_sediments.jl")

if isa(arch, CPU)
    @testset "Doctests" begin
        doctest(OceanBioME)
    end
end