using Documenter, DocumenterCitations, Literate

using OceanBioME
using OceanBioME.SLatissimaModel: SLatissima
using OceanBioME.LOBSTERModel: LOBSTER
using OceanBioME.Boundaries.Sediments: SimpleMultiG, InstantRemineralisation
using OceanBioME.Boundaries: OCMIP_default, GasExchange

using CairoMakie
CairoMakie.activate!(type = "svg")

include("display_parameters.jl")

bib_filepath = joinpath(dirname(@__FILE__), "oceanbiome.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "Simple column model" => "column",
    "Data forced column model" => "data_forced",
    "Model with particles (kelp) interacting with the biogeochemistry" => "kelp",
    "Box model" => "box",
    "Baroclinic instability" => "eady"
]

example_scripts = [ filename * ".jl" for (title, filename) in examples ]

replace_silly_warning(content) = replace(content, r"┌ Warning:.*\s+└ @ JLD2 ~/\.julia/packages/JLD2/.*/reconstructing_datatypes\.jl.*\n" => "")

for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(example_filepath, OUTPUT_DIR; 
                          flavor = Literate.DocumenterFlavor(),
                          repo_root_url = "https://oceanbiome.github.io/OceanBioME.jl",
                          execute = false,
                          postprocess = replace_silly_warning)
    end
end

example_pages = [ title => "generated/$(filename).md" for (title, filename) in examples ]

# create parameter pages

if !isdir(OUTPUT_DIR) mkdir(OUTPUT_DIR) end

model_parameters = (LOBSTER(; grid = BoxModelGrid()),
                    NutrientPhytoplanktonZooplanktonDetritus(; grid = BoxModelGrid()),
                    SLatissima(),
                    TwoBandPhotosyntheticallyActiveRadiation(; grid = BoxModelGrid()),
                    SimpleMultiG(; grid = BoxModelGrid(), depth = 1000),
                    InstantRemineralisation(; grid = BoxModelGrid()),
                    OCMIP_default,
                    GasExchange(; gas = :CO₂).condition.parameters,
                    GasExchange(; gas = :O₂).condition.parameters)

gas_exchange_gas(::Val{G}) where G = G

model_name(model) = if Base.typename(typeof(model)).wrapper == GasExchange
                        "$(gas_exchange_gas(model.gas)) air-sea exchange"
                    else
                        Base.typename(typeof(model)).wrapper
                    end

model_names = [model_name(model) for model in model_parameters]

for (idx, model) in enumerate(model_parameters)
    create_parameter_file!(model, model_names[idx], "$OUTPUT_DIR/$(model_names[idx])_parameters.md")
end

parameter_pages = ["$name" => "generated/$(name)_parameters.md" for name in model_names]

bgc_pages = [
    "Overview" => "model_components/biogeochemical/index.md",
    "LOBSTER" => "model_components/biogeochemical/LOBSTER.md",
    "NPZD" => "model_components/biogeochemical/NPZ.md"
]

sed_pages = [
    "Overview" => "model_components/sediment.md",
]

individuals_pages = [
    "Overview" => "model_components/individuals/index.md",
    "Sugar Kelp (Broch and Slagstad 2012 ++)" => "model_components/individuals/slatissima.md",
]

component_pages = [
    "Biogeochemical models" => bgc_pages,
    "Air-sea gas exchange" => "model_components/air-sea-gas.md",
    "Sediment models" => sed_pages,
    "Light attenuation models" => "model_components/light.md",
    "Individuals" => individuals_pages,
    "Utilities" => "model_components/utils.md"
]

numerical_pages = [
    "Positivity preservation" => "numerical_implementation/positivity-preservation.md"
]

appendix_pages = [
    "Library" => "appendix/library.md",
    "Function index" => "appendix/function_index.md",
    "Parameters" => parameter_pages
]

pages = [
    "Home" => "index.md",
    "Quick start" => "quick_start.md",
    "Model components and setup" => component_pages,
    "Visualization" => "visualization.md",
    "Examples" => example_pages,
    "Numerical implementation" => numerical_pages,
    "References" => "references.md",
    "Appendix" => appendix_pages
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
    collapselevel = 1,
    prettyurls = get(ENV, "CI", nothing) == "true",
    canonical = "https://OceanBioME.github.io/OceanBioME/stable/",
    mathengine = MathJax3(),
    assets = String["assets/citations.css"]
)

makedocs(bib,
    sitename = "OceanBioME.jl",
    authors = "Jago Strong-Wright, John R. Taylor, and Si Chen",
    format = format,
    pages = pages,
    modules = [OceanBioME],
    doctest = false,
    draft = true,
    strict = true,
    clean = true,
    checkdocs = :exports
)

@info "Clean up temporary .jld2/.nc files created by doctests..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

deploydocs(
    repo = "github.com/OceanBioME/OceanBioME.jl",
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
    forcepush = true,
    push_preview = true,
    devbranch = "main"
)
