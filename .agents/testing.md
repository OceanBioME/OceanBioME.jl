# Testing Guidelines

## Layout

- Tests are flat files `test/test_<area>.jl`, each `include`d from `test/runtests.jl`.
  A new test file must be added to `runtests.jl` or it never runs
- `test/dependencies_for_runtests.jl` loads `OceanBioME, Test, CUDA, Oceananigans, JLD2,
  Oceananigans.Units, Documenter` and defines the global `architecture`
- Box-model tests only run on CPU; doctests (`doctest(OceanBioME)`) run at the end, on CPU
- Test-only dependencies are in `[extras]` / `[targets] test` in the root `Project.toml`
- `test/test_sugar_kelp.jl` exists but is not currently included in `runtests.jl`

## Running Tests

`dependencies_for_runtests.jl` errors unless CUDA is functional or `CUDA_VISIBLE_DEVICES=-1`
is set, so on machines without an NVIDIA GPU always set it.

```sh
# Full suite on CPU (there is no test-group selection)
CUDA_VISIBLE_DEVICES=-1 DATADEPS_ALWAYS_ACCEPT=true julia --project -e 'using Pkg; Pkg.test()'
```

A single file needs the test environment active, e.g. with TestEnv.jl installed in the global
environment:

```sh
CUDA_VISIBLE_DEVICES=-1 DATADEPS_ALWAYS_ACCEPT=true julia --project -e '
using TestEnv; TestEnv.activate()
include("test/dependencies_for_runtests.jl")
include("test/test_light.jl")
'
```

`test_solvers.jl` and `test_external_plankton_interface.jl` don't include
`dependencies_for_runtests.jl` themselves, hence including it first.

Doctests alone:

```sh
CUDA_VISIBLE_DEVICES=-1 julia --project -e '
using TestEnv; TestEnv.activate()
using Documenter, OceanBioME
doctest(OceanBioME)
'
```

## Which tests to run

| Changed area | Test files |
|---|---|
| `src/Models/AdvectedPopulations/NutrientsPlanktonDetritus/` | `test_NutrientsPlanktonDetritus.jl`, `test_external_plankton_interface.jl`, `test_construction_in_functions.jl` |
| `src/Models/AdvectedPopulations/PISCES/` | `test_PISCES.jl` |
| `src/Light/` | `test_light.jl` |
| `src/Models/CarbonChemistry/`, `src/Models/GasExchange/` | `test_gasexchange_carbon_chem.jl` |
| `src/Sediments/`, `src/Models/Sediments/` | `test_sediments.jl` |
| `src/Particles/`, `src/Models/Individuals/` | `test_particles.jl`, `test_sugar_kelp.jl` |
| `src/BoxModel/` | `test_boxmodel.jl` |
| `src/Utils/` | `test_utils.jl`, `test_solvers.jl` |
| `src/OceanBioME.jl` | `test_construction_in_functions.jl`, then everything |

## CI

Tests run on Buildkite (`.buildkite/pipeline.yml`) on CPU and on a CUDA GPU, with a single Julia
version. Results appear on PRs as the `CPU tests` / `GPU tests` statuses (`gh pr checks`).

## Writing Tests

- Start the file with `include("dependencies_for_runtests.jl")` and build grids with
  `RectilinearGrid(architecture; ...)` so it runs on CPU and GPU
- Group with nested `@testset`s
- For new biogeochemistry, test conservation of the model's currency (e.g. total nitrogen)
- Test numerical accuracy where analytical solutions exist
- Test `Float32` construction/time-stepping where practical
- Use minimal grid sizes and few time steps; avoid hardcoded indices (use `size(grid, d)`)
- Avoid `@allowscalar` — use `on_architecture(CPU(), interior(field))`

## Debugging Tips

- GPU "dynamic invocation error": run on CPU first. If the error goes away, the problem is
  GPU-specific — usually a missing `@inline`, type instability, or a struct without
  `Adapt.adapt_structure`
- Julia version / resolver issues: delete `Manifest.toml` and run `using Pkg; Pkg.instantiate()`
