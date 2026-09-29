---
name: add-feature
description: Checklist for adding new biogeochemistry, light, sediment, particle, or utility features to OceanBioME
---

# Add Feature

## Step 0: Decide where it goes

- **New plankton / nutrient / detritus / carbon / oxygen behaviour** that fits
  nutrients → plankton → detritus: write an **NPD component** under
  `src/Models/AdvectedPopulations/NutrientsPlanktonDetritus/<Component>/`
  (see `docs/src/model_implementation.md`). Plankton components implement
  `nutrient_uptake`, `dissolved_waste`, `solid_waste`, `inorganic_waste`, and the framework
  closes the budget
- **A structurally different model**: subtype `AbstractContinuousFormBiogeochemistry` or
  `AbstractBiogeochemistry` (see `docs/src/model_implementation_custom.md`, and `PISCES/`)
- **Light**: `src/Light/`; **sediments**: framework in `src/Sediments/`, biogeochemistry in
  `src/Models/Sediments/`; **individuals**: `src/Particles/` + `src/Models/Individuals/`;
  **carbonate chemistry / gas exchange**: `src/Models/CarbonChemistry/`, `src/Models/GasExchange/`
- **Something that modifies the state after tendencies/updates** (like `ScaleNegativeTracers`):
  a `modifiers` entry, in `src/Utils/`

Look at the closest existing sibling and copy its structure.

## Checklist

1. **Get the source**: equations and parameter values from a paper the user supplies. Don't
   invent parameters or references; ask if they're missing
2. **Define the struct**, parametric in `FT`, with a hand-written docstring listing keyword
   arguments **with units**
3. **Constructor(s)**: `Component(FT = Float64; kwargs...)` with literal defaults commented with
   units, converted with `convert(FT, ·)`; add `Component(grid::AbstractGrid{FT}; ...)` if it
   needs fields (e.g. sinking velocities via `setup_velocity_fields`)
4. **Oceananigans interface**: `required_biogeochemical_tracers`,
   `required_biogeochemical_auxiliary_fields`, `biogeochemical_drift_velocity` as needed —
   `import` them, don't just define them
5. **Tendencies**: `@inline`, type-stable, `ifelse` not `if`, `zero(FT)` not `0.0`, dispatch on
   `Val{:tracer}`; reference the source equation in a comment
6. **GPU**: `Adapt.adapt_structure` if the struct holds fields/arrays
7. **`show` / `summary`** methods consistent with neighbours
8. **Exports** in `src/OceanBioME.jl` if user-facing; exported names need docstrings
   (`checkdocs = :exports`); a new submodule needs an `@autodocs` block in
   `docs/src/appendix/library.md`
9. **Tests** in the matching `test/test_*.jl` (or a new file added to `runtests.jl`):
   construction, a few time steps on `architecture`, conservation of the currency, `Float32`
10. **Docs**: a page under `docs/src/model_components/` added to `pages` in `docs/make.jl`;
    add to `model_parameters` in `docs/make.jl` if it should get a parameter page; bibtex in
    `docs/oceanbiome.bib` only from supplied/retrieved sources
11. **Validation** script in `validation/` when reproducing published results
12. **Run** `/run-tests` for the affected files on CPU

## Key Conventions

- File names: snake_case; module files PascalCase matching the module
- Type names: PascalCase; function names: snake_case
- Tracer names: reuse existing ones (`NO₃`, `NH₄`, `P`, `Z`, `DIC`, `Alk`, `O₂`, ...)
- `grid` is positional, keyword args after the semicolon
