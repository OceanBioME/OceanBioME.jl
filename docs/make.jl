pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add OceanBioME to environment stack

using Documenter, DocumenterCitations, Literate, Glob

using OceanBioME
using OceanBioME.SLatissimaModel: SLatissima
using OceanBioME.LOBSTERModel: LOBSTER
using OceanBioME.Boundaries.Sediments: SimpleMultiG

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
    "eady.jl"
]

silly_warnings = ["┌ Warning: type parameters for NamedTuple{(:architecture, :growth_rate_adjustement, :photosynthetic_efficiency, :minimum_carbon_reserve, :structural_carbon, :exudation, :erosion, :saturation_irradiance, :structural_dry_weight_per_area, :structural_dry_to_wet_weight, :carbon_reserve_per_carbon, :nitrogen_reserve_per_nitrogen, :minimum_nitrogen_reserve, :maximum_nitrogen_reserve, :growth_adjustement_2, :growth_adjustement_1, :maximum_specific_growth_rate, :structural_nitrogen, :photosynthesis_at_ref_temp_1, :photosynthesis_at_ref_temp_2, :photosynthesis_ref_temp_1, :photosynthesis_ref_temp_2, :photoperiod_1, :photoperiod_2, :respiration_at_ref_temp_1, :respiration_at_ref_temp_2, :respiration_ref_temp_1, :respiration_ref_temp_2, :photosynthesis_arrhenius_temp, :photosynthesis_low_temp, :photosynthesis_high_temp, :photosynthesis_high_arrhenius_temp, :photosynthesis_low_arrhenius_temp, :respiration_arrhenius_temp, :current_speed_for_0p65_uptake, :nitrate_half_saturation, :ammonia_half_saturation, :maximum_nitrate_uptake, :maximum_ammonia_uptake, :current_1, :current_2, :current_3, :respiration_reference_A, :respiration_reference_B, :exudation_redfield_ratio, :pescribed_velocity, :pescribed_temperature, :pescribed_salinity, :x, :y, :z, :A, :N, :C, :nitrate_uptake, :ammonia_uptake, :primary_production, :frond_exudation, :nitrogen_erosion, :carbon_erosion, :custom_dynamics, :scalefactor, :latitude),Tuple} do not match type NamedTuple in workspace; reconstructing
└ @ JLD2 ~/.julia/packages/JLD2/ryhNR/src/data/reconstructing_datatypes.jl:475",
"┌ Warning: type Oceananigans.TurbulenceClosures.ScalarDiffusivity{Oceananigans.TurbulenceClosures.ExplicitTimeDiscretization,Oceananigans.TurbulenceClosures.ThreeDimensionalFormulation,Main.##420.#κₜ,NamedTuple{(:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM, :DIC, :Alk),Tuple{Vararg{T, N}} where {N, T}{9,Main.##420.#κₜ}}} does not exist in workspace; reconstructing
└ @ JLD2 ~/.julia/packages/JLD2/ryhNR/src/data/reconstructing_datatypes.jl:403",
"┌ Warning: type Oceananigans.TurbulenceClosures.ScalarDiffusivity{Oceananigans.TurbulenceClosures.ExplicitTimeDiscretization,Oceananigans.TurbulenceClosures.ThreeDimensionalFormulation,Main.##318.#κₜ,NamedTuple{(:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM, :DIC, :Alk),Tuple{Vararg{T, N}} where {N, T}{9,Main.##318.#κₜ}}} does not exist in workspace; reconstructing
└ @ JLD2 ~/.julia/packages/JLD2/ryhNR/src/data/reconstructing_datatypes.jl:403"]

function replace_silly_warning(content)
    for str in silly_warinings
        content = replace(content, str => "")
    end
    return content
end

for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)

    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(example_filepath, OUTPUT_DIR; 
                          flavor = Literate.DocumenterFlavor(), 
                          repo_root_url="https://oceanbiome.github.io/OceanBioME.jl", 
                          execute=true,
                          postprocess=replace_silly_warning)
    end
end

example_pages = [
    "Simple column model" => "generated/column.md",
    "Data forced column model" => "generated/data_forced.md",
    "Model with particles (kelp) interacting with the biogeochemistry" => "generated/kelp.md",
    "Box model" => "generated/box.md",
    "Baroclinical instability" => "generated/eady.md"]

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


param_pages = [
    "Overview" => "appendix/params/index.md",
    "LOBSTER" => "appendix/params/LOBSTER.md",
    "SLatissima" => "appendix/params/SLatissima.md"
]

appendix_pages = [
    "Library" => "appendix/library.md",
    "Function index" => "appendix/function_index.md",
    "Parameters" => param_pages
]

pages = [
    "Home" => "index.md",
    "Quick start" => "quick_start.md",
    "Model components and setup" => component_pages,
    "Examples" => example_pages,
    "Numerical implementation" => numerical_pages,
    "Gallery" => "gallery.md",
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

deploydocs(
          repo = "github.com/OceanBioME/OceanBioME.jl",
          forcepush = true,
          push_preview = true,
          devbranch = "main"
)

