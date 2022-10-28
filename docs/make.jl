using Documenter, DocumenterCitations, Literate, Glob
using OceanBioME

bib_filepath = joinpath(dirname(@__FILE__), "oceanbiome.bib")
bib = CitationBibliography(bib_filepath)

# Examples

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

examples = [
    "box.jl",
    "column.jl",
    "data_forced.jl",
    "kelp.jl",
    "sediment.jl"
]
# going to need some work done on the formatting of the comments in examples to get them to render properly
#This *runs* all of the example so can't really use this (currently have #example replaced with #example so they won't run)
#=for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end
=#
example_pages = [
    "Simple virtual aquarium (column model)" => "generated/column.md",
    "Data forced virtual aquarium" => "generated/data_forced.md",
    "Virtual aquarium with sediment" => "generated/sediment.md",
    "Virtual aquarium with kelp" => "generated/kelp.md",
    "Box model" => "generated/box.md"
]

bgc_pages = [
    "Overview" => "model_components/biogeochemical/index.md",
    "PISCES" => "model_components/biogeochemical/PISCES.md",
    "LOBSTER" => "model_components/biogeochemical/LOBSTER.md",
    "NPZ" => "model_components/biogeochemical/NPZ.md"
]

sed_pages = [
    "Overview" => "model_components/sediment/index.md",
    "Wang et al. 2020" => "model_components/sediment/wang.md",
    "Soetaert et al. 2000" => "model_components/sediment/soetaert.md",
]

stat_pop_pages = [
    "Overview" => "model_components/quasi-stationary-populations/index.md",
    #"Seagrass (X et al. Y)" => "model_components/quasi-stationary-populations/seagrass.md",
    #"Nektons (Z et al. A)" => "model_components/quasi-stationary-populations/nekton.md",
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
    "Quasi-stationary populations" => stat_pop_pages,
    "Individuals" => individuals_pages,
    "Utilities" => "model_components/utils.md"
]

numerical_pages = [
    "Individuals" => "numerical_implimentation/individuals.md",
    "Sediments" => "numerical_implimentation/sediments.md",
    "Forced auxiliary fields" => "numerical_implimentation/forced-aux-fields.md"
]

setup_pages = [
    "Overview" => "setup/index.md",
    #...
]

param_pages = [
    "Overview" => "appendix/params/index.md",
    "PISCES" => "appendix/params/PISCES.md",
    "LOBSTER" => "appendix/params/LOBSTER.md",
    "SLatissima" => "appendix/params/SLatissima.md"
]

appendix_pages = [
    "Function index" => "appendix/function_index.md",
    "Parameters" => param_pages
]

pages = [
    "Home" => "index.md",
    "Quick start" => "quick_start.md",
    "Examples" => example_pages,
    "Model components" => component_pages,
    "Numerical implimentation" => numerical_pages,
    #"Model setup" => setup_pages,
    "Notes" => "notes.md",
    "Gallery" => "gallery.md",
    "References" => "references.md",
    "Appendix" => appendix_pages
]

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
    collapselevel = 1,
    prettyurls = false,#get(ENV, "CI", nothing) == "true",
    canonical = "https://OceanBioME.github.io/OceanBioME/stable/",
    mathengine = MathJax3(),
)

makedocs(bib,
    sitename = "OceanBioME.jl",
    authors = "Jago Strong-Wrigt, John R. Taylor, and Si Chen",
    format = format,
    pages = pages,
    modules = Module[OceanBioME],
    doctest = false,#true,
    strict = false,#true,
    clean = true,
    checkdocs = :exports
)

@info "Cleaning up temporary .jld2 and .nc files created by doctests..."

for file in vcat(glob("docs/*.jld2"), glob("docs/*.nc"))
    rm(file)
end
#=
deploydocs(
          repo = "github.com/OceanBioME/OceanBioME/Documentation.git",
      versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
     forcepush = true,
  push_preview = true,
     devbranch = "main"
)
=#
