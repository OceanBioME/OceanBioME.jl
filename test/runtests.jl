include("dependencies_for_runtests.jl")

group = get(ENV, "TEST_GROUP", :all) |> Symbol

@info "Please do not merge this branch I just want to see if it breaks the test"

if group == :all
    include("test_utils.jl")
    include("test_light.jl")
    include("test_slatissima.jl")
    include("test_LOBSTER.jl")
    include("test_NPZD.jl")
    include("test_gasexchange.jl")
    include("test_sediments.jl")

    if architecture == CPU() 
        include("test_boxmodel.jl") # box models (probably) don't work on GPU, and it wouldn't be faster anyway
    end

    architecture == CPU() && @testset "Doctests" begin
        doctest(OceanBioME)
    end
end
