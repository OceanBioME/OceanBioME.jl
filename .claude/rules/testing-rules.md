---
paths:
  - test/**/*.jl
---

# Testing Rules

## Layout

- Tests are flat files `test/test_<area>.jl`, each `include`d from `test/runtests.jl`.
  A new test file must be added to `runtests.jl` or it never runs
- `test/dependencies_for_runtests.jl` loads `OceanBioME, Test, CUDA, Oceananigans, JLD2,
  Oceananigans.Units, Documenter` and defines the global `architecture`
- Box-model tests only run on CPU; doctests (`doctest(OceanBioME)`) run at the end on CPU
- `test_sugar_kelp.jl` exists but is not currently included in `runtests.jl`
- Test-only dependencies are in `[extras]` / `[targets] test` in the root `Project.toml`

## Architecture selection

`dependencies_for_runtests.jl` **errors** unless CUDA is functional or
`CUDA_VISIBLE_DEVICES=-1` is set. On a machine without an NVIDIA GPU (e.g. macOS), always set
`CUDA_VISIBLE_DEVICES=-1`.

## Running Tests

There is no test-group selection: `Pkg.test()` runs everything.

```sh
# Full suite on CPU
CUDA_VISIBLE_DEVICES=-1 julia --project -e 'using Pkg; Pkg.test()'
```

To run a single file you need the test environment active. With
[TestEnv.jl](https://github.com/JuliaTesting/TestEnv.jl) installed in the global environment:

```sh
CUDA_VISIBLE_DEVICES=-1 julia --project -e '
using TestEnv; TestEnv.activate()
include("test/dependencies_for_runtests.jl")
include("test/test_light.jl")
'
```

(Including `dependencies_for_runtests.jl` first is needed for `test_solvers.jl` and
`test_external_plankton_interface.jl`, which don't include it themselves; including it twice is harmless.)

## Writing Tests

- Start the file with `include("dependencies_for_runtests.jl")` and build grids with
  `RectilinearGrid(architecture; ...)` so the file runs on CPU and GPU
- Name test files `test_<area>.jl`; group with nested `@testset`s
- Test conservation of the model's currency (e.g. total nitrogen/carbon) for new biogeochemistry
- Test numerical accuracy where analytical solutions exist (e.g. PAR attenuation, box model decay)
- Test `Float32` construction/time-stepping for new components where practical
- Use minimal grid sizes and few time steps to keep CI time down
- Avoid hardcoded grid indices — use `size(grid, d)` instead of literal numbers

## Quality

- **Avoid `@allowscalar` in new tests** — move data to CPU with
  `on_architecture(CPU(), interior(field))` first
- Always add tests for new functionality
- Doctests in `src/` are run by `Pkg.test()` on CPU, so a changed `show` method can fail the suite

## Debugging

- GPU "dynamic invocation error": run on CPU first to isolate GPU-specific issues, then look for
  a missing `@inline`, a type instability, or a struct missing `adapt_structure`
- Julia version / resolver issues: delete `Manifest.toml`, then `Pkg.instantiate()`
- Test data is fetched with DataDeps (`datadep"test_data/..."`); set `DATADEPS_ALWAYS_ACCEPT=true`
  for non-interactive runs
