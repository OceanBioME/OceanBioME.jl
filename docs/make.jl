using Documenter, DocumenterCitations, Literate

using OceanBioME
using OceanBioME: SugarKelp, LOBSTER, NutrientPhytoplanktonZooplanktonDetritus, SimpleMultiGSediment, InstantRemineralisationSediment
using OceanBioME: CarbonChemistry, GasExchange

using Oceananigans.Grids: RectilinearGrid

using CairoMakie
CairoMakie.activate!(type = "svg")

include("display_parameters.jl")

bib_filepath = joinpath(dirname(@__FILE__), "oceanbiome.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "Box model" => "box",
    "Simple column model" => "column",
    "Baroclinic instability" => "eady",
    "Data forced column model" => "data_forced",
    "Model with particles (kelp) interacting with the biogeochemistry" => "kelp",
    "Data assimilation" => "data_assimilation"
]

example_scripts = [ filename * ".jl" for (title, filename) in examples ]

replace_silly_warning(content) = replace(content, r"┌ Warning:.*\s+└ @ JLD2 ~/\.julia/packages/JLD2/.*/reconstructing_datatypes\.jl.*\n" => "")

Threads.@threads for example in example_scripts
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(example_filepath, OUTPUT_DIR; 
                          flavor = Literate.DocumenterFlavor(),
                          repo_root_url = "https://oceanbiome.github.io/OceanBioME.jl",
                          execute = true,
                          postprocess = replace_silly_warning)
    end
end

example_pages = [ title => "generated/$(filename).md" for (title, filename) in examples ]

# create parameter pages

if !isdir(OUTPUT_DIR) mkdir(OUTPUT_DIR) end

small_grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))

model_parameters = (LOBSTER(; grid = BoxModelGrid(), light_attenuation_model = nothing).underlying_biogeochemistry,
                    NutrientPhytoplanktonZooplanktonDetritus(; grid = BoxModelGrid(), light_attenuation_model = nothing).underlying_biogeochemistry,
                    SugarKelp(),
                    TwoBandPhotosyntheticallyActiveRadiation(; grid = small_grid),
                    SimpleMultiGSediment(small_grid).biogeochemistry,
                    InstantRemineralisationSediment(small_grid).biogeochemistry,
                    CarbonChemistry(),
                    CarbonDioxideGasExchangeBoundaryCondition().condition.func,
                    OxygenGasExchangeBoundaryCondition().condition.func)

exchanged_gas(::Val{G}) where G = G

model_name(model) = if Base.typename(typeof(model)).wrapper == GasExchange
                        ifelse(isa(model.water_concentration, CarbonChemistry), "CO₂", "O₂")*" air-sea exchange"
                    else
                        Base.typename(typeof(model)).wrapper
                    end

model_names = [model_name(model) for model in model_parameters]

for (idx, model) in enumerate(model_parameters)
    create_parameter_file!(model, model_names[idx], "$OUTPUT_DIR/$(model_names[idx])_parameters.md")
end

parameter_pages = ["$name" => "generated/$(name)_parameters.md" for name in model_names]

pisces_pages = ["PISCES" => "model_components/biogeochemical/PISCES/PISCES.md",
                "Queries" => "model_components/biogeochemical/PISCES/notable_differences.md"]

bgc_pages = [
    "Overview" => "model_components/biogeochemical/index.md",
    "LOBSTER" => "model_components/biogeochemical/LOBSTER.md",
    "NPZD" => "model_components/biogeochemical/NPZ.md",
    "PISCES" => pisces_pages
]

sediments_pages = [
    "Overview" => "model_components/sediments/index.md",
    "Simple Multi-G" => "model_components/sediments/simple_multi_g.md",
    "Instant remineralisation" => "model_components/sediments/instant_remineralisation.md"
]

individuals_pages = [
    "Overview" => "model_components/individuals/index.md",
    "Sugar Kelp (Broch and Slagstad 2012 ++)" => "model_components/individuals/slatissima.md",
]

component_pages = [
    "Biogeochemical models" => bgc_pages,
    "Air-sea gas exchange" => "model_components/air-sea-gas.md",
    "Carbon chemistry" => "model_components/carbon-chemistry.md",
    "Sediment models" => sediments_pages,
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
    "Examples" => example_pages,
    "Model components and setup" => component_pages,
    "Implementing new models" => "model_implementation.md",
    "Numerical implementation" => numerical_pages,
    "Visualization" => "visualization.md",
    "Contibutors guide" => "contributing.md",
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

makedocs(sitename = "OceanBioME.jl",
         authors = "Jago Strong-Wright and contributors",
         format = format,
         pages = pages,
         modules = [OceanBioME],
         plugins = [bib],
         doctest = false,#true,
         clean = true,
         checkdocs = :exports)

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
    versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
    forcepush = true,
    push_preview = true,
    devbranch = "main"
)
