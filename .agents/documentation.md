# Documentation

## Building Docs Locally

Examples are run separately (as on CI) and write markdown into `docs/src/generated/`, which
`docs/make.jl` then uses. `docs/make.jl` also generates the model parameter pages.

```sh
# Set up the docs environment
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

# Run one example (name without .jl) or all of them
julia --project=docs/ docs/make_examples.jl column
julia docs/prepare_examples.jl

# Build
JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
```

Then open `docs/build/index.html`, or `using LiveServer; serve(dir="docs/build")`.

## Fast Local Builds

Running all examples is slow (`eady`, `kelp` and `oae_experiment` are the heavy ones). To check
prose and cross-references only, temporarily edit `docs/make.jl`:

1. Remove the unbuilt entries from `examples` (pages point at `generated/<name>.md`, which only
   exists once the example has run)
2. If needed, add `warnonly = [:cross_references, :missing_docs]` to `makedocs`
3. Optionally `draft = true`

**Revert these changes before committing.** Note `doctest = false` is already set in
`docs/make.jl`; doctests are run by `Pkg.test()` instead.

## Documentation Style

- Documenter.jl cross-references: `[text](@ref id)` with `@id` anchors on headers
- Citations via DocumenterCitations: add bibtex to `docs/oceanbiome.bib` and cite with
  `(@citet)` / `(@citep)`. Only add references supplied by a user or actually retrieved; never
  write bibtex from memory
- LaTeX (```` ```math ````) is fine in `docs/src` pages; state units in parameter tables
- New pages must be added to `pages` in `docs/make.jl`
- Exported names need docstrings (`checkdocs = :exports`); a new submodule needs an `@autodocs`
  block in `docs/src/appendix/library.md`
- New models with parameters can be added to `model_parameters` in `docs/make.jl` to get a
  generated parameter page

## Docstrings

OceanBioME does not use DocStringExtensions; write signatures by hand:

~~~~
"""
    PhytoZoo([FT = Float64;] kwargs...)

Description.

Keyword Arguments
=================

- `maximum_grazing_rate`: maximum zooplankton grazing rate (1/s)
"""
~~~~

Prefer `jldoctest` over plain `julia` blocks — doctests are run by `Pkg.test()` on CPU, plain
blocks rot. Existing example (`src/Models/CarbonChemistry/carbon_chemistry.jl`):

~~~~
```jldoctest
julia> using OceanBioME

julia> carbon_chemistry = CarbonChemistry()
`CarbonChemistry` model which solves for pCO₂ and pH

julia> pCO₂ = carbon_chemistry(; DIC = 2000.0, Alk = 2000.0, T = 10.0, S = 35.0)
1308.1006995915372
```
~~~~

- Always include the expected output, copied from an actual run
- Exercise `Base.show` where possible; don't end on a boolean comparison
- Use `filter = r"..."` for varying output (see `src/Models/Sediments/simple_multi_G.jl`)
- Use unicode (`μ`, `Δt`) rather than LaTeX in docstrings

## Writing Examples

- Start with `# # [Title](@id anchor)`, explain what the simulation does and which features it
  demonstrates, and include the commented "Install dependencies" block
- Literate style: let code speak for itself. Single `#` → markdown, `##` → code comment
- `using OceanBioME, Oceananigans` (+ `Oceananigans.Units`); only import unexported names
  explicitly
- State units of forcing and parameters; use `nothing #hide` to suppress clutter
- Only packages in `docs/Project.toml` are available on CI
- To publish an example it must be listed in the Buildkite matrix in `.buildkite/pipeline.yml`,
  in `examples` in `docs/make.jl`, and in `EXAMPLES` in `docs/prepare_examples.jl`
