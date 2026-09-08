using Pkg

# Setup OceanBioME as a development dependency
# @__DIR__ is docs/, so ".." goes to OceanBioME.jl
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))

const EXAMPLES = [
    "box",
    "column",
    "eady",
    "kelp",
    "data_assimilation",
    "oae_experiment"
]

for example_name in EXAMPLES
    @info "Processing example: $example_name"
    
    try
        run(`julia --project=docs/ docs/make_examples.jl $example_name`)
    catch e
        @error "Failed to process: $example_name" exception=(e, catch_backtrace())
    end
end

@info "All examples processed!"
