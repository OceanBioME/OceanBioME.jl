---
paths:
  - examples/**/*.jl
---

# Examples Rules

## Writing Examples

- Start with a `# # [Title](@id anchor)` header and explain what the simulation does and which
  OceanBioME features it demonstrates
- Include the "Install dependencies" block (commented `pkg"add ..."`) like the existing examples
- Let code "speak for itself" — keep explanations concise (Literate style)
- New examples should add value while remaining simple, and should run in reasonable time on CI
- Load with `using OceanBioME, Oceananigans` (plus `Oceananigans.Units`); only import
  unexported names explicitly (`using Oceananigans.Fields: FunctionField`)
- State units of forcing functions and parameters in comments
- End setup-only lines with `nothing #hide` where output would clutter the page
- Only packages in the docs environment (`docs/Project.toml`) are available on CI

## Adding an example to the docs

An example is only built if it is listed in **all** of:
- the `matrix.setup.example_path` list in `.buildkite/pipeline.yml`
- `examples` in `docs/make.jl`
- `EXAMPLES` in `docs/prepare_examples.jl`

(`data_forced.jl` is currently not in these lists.)

## Literate.jl Comment Conventions

- Single `#` comments become markdown blocks in generated documentation
- Double `##` comments remain as code comments within code blocks
- Use single `#` only for narrative text that should render as markdown
