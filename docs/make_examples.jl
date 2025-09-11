using Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

replace_silly_warning(content) = replace(content, r"┌ Warning:.*\s+└ @ JLD2 ~/\..*/packages/JLD2/.*/reconstructing_datatypes\.jl.*\n" => "")

example_name = ARGS[1]

example_filepath = joinpath(EXAMPLES_DIR, example_name*".jl")

Literate.markdown(example_filepath, OUTPUT_DIR; 
                  flavor = Literate.DocumenterFlavor(),
                  repo_root_url = "https://oceanbiome.github.io/OceanBioME.jl",
                  execute = true,
                  postprocess = replace_silly_warning)
