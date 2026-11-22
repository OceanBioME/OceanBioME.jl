# OceanBioME.jl — Agent Rules

## Project Overview

OceanBioME.jl is a Julia package for modelling the coupled interactions between ocean
biogeochemistry, carbonate chemistry, and physics. It is built on top of
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and plugs into Oceananigans models
through the `biogeochemistry` keyword (`Oceananigans.Biogeochemistry` interface), and also provides
a standalone `BoxModel` for 0D integration.

Main components:
- **Biogeochemical models**: the `NutrientsPlanktonDetritus` (NPD) framework with composable
  `nutrients` / `plankton` / `detritus` / `inorganic_carbon` / `oxygen` components, and presets
  (`NPZD`, `LOBSTER`, `ImplicitBiology`); `PISCES`
- **Light attenuation**: two-band, multi-band, prescribed PAR
- **Carbon chemistry and air–sea gas exchange**: `CarbonChemistry`, `GasExchange`
- **Sediments**: `SimpleMultiGSediment`, `InstantRemineralisationSediment`
- **Individuals / particles**: `SugarKelp`, `BiogeochemicalParticles`
- **Utilities**: positivity preservation (`ScaleNegativeTracers`, `ZeroNegativeTracers`), `Budget`,
  sinking-velocity fields, timestep helpers

## Language & Environment

- **Julia 1.10+** (CI runs a newer Julia, see `.buildkite/pipeline.yml`) | CPU and CUDA GPU
- **Key packages**: Oceananigans.jl, KernelAbstractions.jl, CUDA.jl, Adapt.jl
- Only CUDA GPUs are tested. Do not add Metal/AMD/Reactant/Enzyme-specific code unless asked

## Critical Rules

### Kernel and tendency functions (GPU compatibility)

Most OceanBioME physics runs **inside Oceananigans' tracer tendency kernels**: biogeochemical
tendencies are functions called per grid point, e.g.
`(bgc::NPD)(i, j, k, grid, val_name::Val, clock, fields, auxiliary_fields)` or, for
continuous-form models, `(bgc)(::Val{:P}, x, y, z, t, tracers...)`. Treat these as kernel code:

- Mark them and everything they call `@inline`
- Type-stable, allocation-free, no error messages, no `Model`s
- Use `ifelse`, not short-circuiting `if`/`else`, on values that vary per grid point
- Select tracer behaviour by dispatch on `Val{:name}`, not by comparing symbols at runtime
- Where OceanBioME launches its own kernels (sediments, particles, negative tracers), use
  `@kernel` / `@index` and `launch!` from `Oceananigans.Utils`
- **Never loop over grid points outside kernels**
- Structs used in kernels must have an `Adapt.adapt_structure` method if they hold arrays or fields

### Type Stability & Floating-Point Type

- All structs must be concretely typed and parametric in `FT` where they hold parameters
- Type annotations are for **dispatch**, not documentation
- **Never hardcode Float64** in kernels: use `zero(FT)`, `one(FT)`, `convert(FT, 1//2)`, or
  rational literals
- Constructors follow the pattern `Component(FT = Float64; kwargs...)` and/or
  `Component(grid::AbstractGrid{FT}; kwargs...)`; literal keyword defaults are fine and are
  converted with `convert(FT, value)` when building the struct

### Parameters and units

- Every dimensional parameter default gets a trailing comment with its units, e.g.
  `maximum_grazing_rate = 9.26e-6, # 1/s`. Rates are per second (SI time); concentrations are
  usually mmol/m³ of the model's currency (e.g. mmol N/m³)
- State the source of parameter values (paper, fitted to what) in a comment or docstring when known
- Do not invent parameter values or literature references. If a value's source is unknown, say so

### Imports

- Source code: explicit imports (`using Oceananigans.Fields: ZeroField`), `import` for methods
  being extended (`import Oceananigans.Biogeochemistry: required_biogeochemical_tracers`).
  There is no automated ExplicitImports check here, so be careful
- Examples/docs: `using OceanBioME, Oceananigans`; import unexported names explicitly only when
  needed (e.g. `using Oceananigans.Fields: FunctionField`)

### Docstrings

- Write the signature by hand at the top of the docstring (DocStringExtensions is **not** a
  dependency), followed by a description and a `Keyword Arguments` section
- **Prefer `jldoctest` blocks over plain `julia` blocks** — doctests run in the CPU test suite
  via `doctest(OceanBioME)`; plain blocks rot
- Include `# output` with verifiable output; prefer `show` methods over boolean comparisons
- Use unicode for math (`Δt`, `μ`, `PAR`), not LaTeX, in docstrings. LaTeX (```` ```math ````)
  is fine in the `docs/src` markdown pages

### Model Constructors

- OceanBioME biogeochemistry takes `grid` positionally: `NPZD(grid; light_attenuation = ...)`,
  `LOBSTER(grid)`, `PISCES(grid)`
- Oceananigans models: `NonhydrostaticModel(grid; biogeochemistry)`,
  `HydrostaticFreeSurfaceModel(grid; biogeochemistry)`; box models: `BoxModel(; biogeochemistry, clock)`
- Omit the semicolon when there are no keyword arguments: `LOBSTER(grid)` not `LOBSTER(grid;)`

## Naming Conventions

- **Files**: snake_case, matching the type they define where there is one; module files are
  PascalCase and match the module (`Light/Light.jl`, `Sediments/Sediments.jl`)
- **Types/Constructors**: PascalCase **only for true constructors**
- **Functions**: snake_case; functions that return values are never PascalCase
- **Kernels**: may prefix with underscore — `_scale_negative_tracers!`
- **Tracer names**: follow existing conventions (`NO₃`, `NH₄`, `P`, `Z`, `sPOM`, `bPOM`, `DIC`,
  `Alk`, `O₂`); check what the relevant component already uses before inventing a new one
- **Variables**: verbose English in user-facing APIs and keyword arguments
  (`phytoplankton_mortality_rate`), readable unicode math inside tendency functions (`μ`, `Kₙ`)
- Spelling follows the existing code (British: `remineralisation`, `parameterisation`); don't
  rename existing API to "fix" spelling

## Module Structure

```
src/
├── OceanBioME.jl                 # Main module, exports, Biogeochemistry wrapper
├── Light/                        # PAR / light attenuation models
├── Particles/                    # BiogeochemicalParticles and tracer coupling
├── Sediments/                    # Sediment framework (FlatSediment, bottom coupling, timesteppers)
├── BoxModel/                     # 0D BoxModel, timesteppers, SpeedyOutput
├── Models/
│   ├── AdvectedPopulations/
│   │   ├── NutrientsPlanktonDetritus/   # NPD framework and its components
│   │   └── PISCES/
│   ├── CarbonChemistry/
│   ├── GasExchange/
│   ├── Sediments/                # Sediment biogeochemistry (SimpleMultiG, InstantRemineralisation)
│   └── Individuals/SugarKelp/
└── Utils/                        # Negative tracers, Budget, sinking velocities, timestep helpers
test/        # flat test_*.jl files, included from runtests.jl
docs/        # Documenter + Literate; parameter pages generated in docs/make.jl
examples/    # Literate examples (built on CI)
validation/  # validation scripts (not run on CI)
benchmark/
```

## Common Pitfalls

1. **Type instability** in tendency functions ruins GPU performance
2. **Overconstraining types**: use annotations for dispatch, not documentation
3. **Missing `@inline`** on functions called from tendencies — GPU compilation fails or slows
4. **Missing `Adapt.adapt_structure`** for a new struct holding fields/arrays — GPU runs fail
5. **Tracer conservation**: in the NPD framework, write components through the
   `nutrient_uptake` / `dissolved_waste` / `solid_waste` / `inorganic_waste` interface so the
   framework balances the currency. Hand-written tendencies must conserve mass themselves; check it
6. **Subtle bugs from missing method imports**: extending an Oceananigans function without
   `import` silently creates a new, unused function
7. **Extending `getproperty` to fix undefined property bugs**: fix on the caller side instead
8. **"Type is not callable" errors**: variable name shadows a function — rename or qualify
9. **Quick fixes that break correctness**: if a test fails after a change, revisit the original edit
10. **Commented-out code**: delete it. Git is the journal
11. **2D indexing on fields**: always use 3D indexing (`field[i, j, k]`)
12. **Hardcoded Float64** in kernels: use `zero(FT)` etc.
13. **Scope creep in PRs**: keep changes focused on a single concern
14. **Modifying Project.toml dependencies**: never add, remove, or change `[deps]`, `[extras]`, or
    `[targets]` unless the task absolutely requires it. Only touch `[compat]` when explicitly asked.
    The Oceananigans compat bound in particular is tightly pinned and changes need care
15. **Unsolicited validation**: don't add argument checks, warnings, or guards that weren't asked
    for — propose them instead

## Git Workflow

Follow [ColPrac](https://github.com/SciML/ColPrac). Feature branches, descriptive commits,
update tests and docs with code changes, check CI before merging. See `docs/src/contributing.md`.

## Design Principles

- **Dispatch over conditionals**: use Julia's type system and multiple dispatch instead of
  `if`/`else` branching (e.g. `Val{:tracer}` methods, component types in the NPD framework)
- **Compose, don't duplicate**: prefer a new NPD component over a new monolithic model when the
  physics fits the nutrients–plankton–detritus structure
- **Use `on_architecture` for data transfers** — never manual `Array()` / `CuArray()` calls in `src/`
- **Defaults serve the common case**: avoid `nothing` defaults when a concrete default covers
  most usage
- **Keyword argument names must be consistent** across related components and constructors
- **Always use explicit `return`** in functions longer than one expression
- **One operation per line** as default; break long expressions across lines

## Agent Behavior

- Prioritize type stability and GPU compatibility
- Follow established patterns in existing code (look at a sibling component first)
- Add tests for new functionality; update exports in `src/OceanBioME.jl` when adding public API,
  and give exported names docstrings (`docs/make.jl` uses `checkdocs = :exports`). A new
  submodule needs its own `@autodocs` block in `docs/src/appendix/library.md`
- Reference the source equations (paper and equation number) in comments when implementing
  biogeochemistry — only references you actually have, never from memory

## Further Reading

Detailed reference docs are in `.agents/` — read on demand:

| Document | Content |
|----------|---------|
| `.agents/testing.md` | Running and writing tests, which tests to run, debugging |
| `.agents/documentation.md` | Building docs and examples, fast builds, docstrings, writing examples |
| `.agents/validation.md` | Reproducing published results step-by-step |

## Auto-loading Rules (Claude Code)

Rules in `.claude/rules/` load automatically when you touch matching files:
- `kernel-rules.md` — GPU kernel / tendency requirements (`src/`)
- `docstring-rules.md` — docstring and jldoctest conventions (`src/`)
- `style-rules.md` — naming and comment style (`src/`, `test/`, `validation/`, `examples/`)
- `testing-rules.md` — test writing and running (`test/`)
- `docs-rules.md` — documentation building and style (`docs/`)
- `examples-rules.md` — Literate.jl example conventions (`examples/`)
- `julia-repl-rules.md` — prefer an MCP Julia REPL when available

## Skills (Claude Code slash commands)

- `/run-tests` — run targeted tests, prioritized by what's likely to break
- `/build-docs` — build documentation locally
- `/add-feature` — checklist for adding new biogeochemistry, light, sediment, or utility features
- `/new-simulation` — set up, run, and visualize a new OceanBioME simulation
- `/babysit-ci` — monitor Buildkite CI via GitHub statuses, fix small issues, pause on bigger ones
