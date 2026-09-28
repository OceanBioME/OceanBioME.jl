---
paths:
  - docs/**/*
---

# Documentation Rules

## Building Docs

Examples are built separately from the main docs (as on CI), and the generated markdown in
`docs/src/generated/` is then picked up by `docs/make.jl`:

```sh
# One-off: set up the docs environment
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

# Build a single example (name without .jl), or all of them
julia --project=docs/ docs/make_examples.jl column
julia docs/prepare_examples.jl

# Build the docs
JULIA_DEBUG=Documenter julia --project=docs/ docs/make.jl
```

`docs/make.jl` also generates the parameter pages (`docs/src/generated/*_parameters.md`) from
`display_parameters.jl`, so a new model with parameters listed there must construct on a
`BoxModelGrid()` / small `RectilinearGrid`.

## Fast Local Builds

Examples take a long time (some want a GPU). For a prose/cross-reference check, build without
running the examples. Pages point at `generated/<example>.md`, which only exists after the example
has been built, so temporarily:
1. Remove the entries you haven't built from `examples` in `docs/make.jl`
2. If needed, add `warnonly = [:cross_references, :missing_docs]` to `makedocs`

Optional: `draft = true`. **Revert these changes before committing!**

## Viewing Docs

Open `docs/build/index.html`, or

```julia
using LiveServer
serve(dir="docs/build")
```

## Style

- Use Documenter.jl syntax for cross-references (`[text](@ref id)`, `@id` anchors)
- Add paper references to `docs/oceanbiome.bib` and cite with DocumenterCitations
  (`[Author2020](@citet)` / `(@citep)`). Only add references the user supplies or that you
  have retrieved — never write bibtex entries from memory
- LaTeX math (```` ```math ````) is fine in `docs/src` pages; state units in parameter tables
- New pages must be added to `pages` in `docs/make.jl`
- New exported names need docstrings (`checkdocs = :exports`) and new submodules need an
  `@autodocs` block in `docs/src/appendix/library.md`
- In example code, don't explicitly import names already exported by `using OceanBioME, Oceananigans`

## Docstrings

- See `.claude/rules/docstring-rules.md`
