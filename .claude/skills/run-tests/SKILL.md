---
name: run-tests
description: Run targeted OceanBioME tests, prioritized by what's likely to break
---

# Run Tests

Run the test files most likely to be affected by recent changes, one at a time, before
(optionally) the full suite. Fix each failure before moving on.

## Step 1: Identify What Changed

Look at staged, unstaged, and recent commits and pick the affected files:

| Changed area | Test files to run first |
|---|---|
| `src/Models/AdvectedPopulations/NutrientsPlanktonDetritus/` | `test_NutrientsPlanktonDetritus.jl`, `test_external_plankton_interface.jl`, `test_construction_in_functions.jl` |
| `src/Models/AdvectedPopulations/PISCES/` | `test_PISCES.jl` |
| `src/Light/` | `test_light.jl` |
| `src/Models/CarbonChemistry/`, `src/Models/GasExchange/` | `test_gasexchange_carbon_chem.jl` |
| `src/Sediments/`, `src/Models/Sediments/` | `test_sediments.jl` |
| `src/Particles/`, `src/Models/Individuals/` | `test_particles.jl` (and `test_sugar_kelp.jl`, not in `runtests.jl`) |
| `src/BoxModel/` | `test_boxmodel.jl` |
| `src/Utils/` | `test_utils.jl`, `test_solvers.jl` (solvers) |
| `src/OceanBioME.jl` (Biogeochemistry wrapper, exports) | `test_construction_in_functions.jl`, then everything |
| Docstrings / `show` methods | doctests (run at the end of `Pkg.test()` on CPU) |

## Step 2: Run the Most Likely File

Always set `CUDA_VISIBLE_DEVICES=-1` unless a working NVIDIA GPU is present —
`test/dependencies_for_runtests.jl` errors otherwise.

A single file needs the test environment. This uses TestEnv.jl from the global environment;
if `using TestEnv` fails, ask the user before installing it (`julia -e 'using Pkg; Pkg.add("TestEnv")'`
modifies their global environment).

```sh
CUDA_VISIBLE_DEVICES=-1 DATADEPS_ALWAYS_ACCEPT=true julia --project -e '
using TestEnv; TestEnv.activate()
include("test/dependencies_for_runtests.jl")
include("test/test_light.jl")
'
```

First runs precompile Oceananigans and can take many minutes: run in the background or with a
long timeout.

## Step 3: Fix and Iterate

1. If a test fails, fix the issue
2. Re-run the same file to confirm the fix
3. Move on to the next most likely file
4. If `show` methods or docstrings changed, run doctests:
   ```sh
   CUDA_VISIBLE_DEVICES=-1 julia --project -e '
   using TestEnv; TestEnv.activate()
   using Documenter, OceanBioME
   doctest(OceanBioME)
   '
   ```

## Step 4: Full Suite (optional)

```sh
CUDA_VISIBLE_DEVICES=-1 DATADEPS_ALWAYS_ACCEPT=true julia --project -e 'using Pkg; Pkg.test()'
```

There is no group selection; this runs every file in `runtests.jl`. GPU tests only run on CI
(Buildkite) unless the user has a CUDA machine.

## Notes

- GPU-only failures ("dynamic invocation", "unsupported call"): reproduce on CPU first; common
  causes are a missing `@inline`, type instability, or a struct lacking `adapt_structure`
- If Julia version issues arise, delete `Manifest.toml` and run `Pkg.instantiate()`
- Report failures with the relevant output; don't claim tests pass without running them
