---
paths:
  - src/**/*.jl
---

# Kernel and Tendency Function Rules

OceanBioME code runs on CPU and CUDA GPUs. Two kinds of code execute on the device:

1. **Biogeochemical tendency functions**, called per grid point from Oceananigans' tracer tendency
   kernels, e.g. `(bgc::NPD)(i, j, k, grid, val_name::Val, clock, fields, auxiliary_fields)`,
   `component_tendency(i, j, k, grid, component, val_name, bgc, clock, fields, auxiliary_fields)`,
   NPD component functions (`nutrient_uptake`, `dissolved_waste`, ...), and continuous-form
   `(bgc)(::Val{:P}, x, y, z, t, tracers...)`
2. **OceanBioME's own kernels** in `Sediments/`, `Particles/`, `Utils/negative_tracers.jl`, etc.

Both must follow the rules below.

## Requirements

- Mark every function called from a tendency or kernel `@inline`
- Keep them **type-stable** and **allocation-free**
- Use `ifelse` instead of short-circuiting `if`/`else` on values that vary per grid point
- No error messages, `@warn`, or `println` inside kernels/tendencies
- Models **never** go inside kernels
- Dispatch on `Val{:tracer_name}` to select per-tracer behaviour; never compare symbols at runtime
- Own kernels use KernelAbstractions.jl (`@kernel`, `@index`) launched with
  `launch!(arch, grid, :xyz | :xy, kernel!, args...)` from `Oceananigans.Utils`
- **Never loop over grid points outside kernels**

## GPU adaptation

- Any struct passed into a kernel that holds fields, arrays, or other adaptable members needs an
  `Adapt.adapt_structure` method (see existing `adapts.jl` / `adapt_show_methods.jl` files)
- Structs that hold only `FT` parameters and singletons adapt automatically

## Type Stability

- All structs must be concretely typed (parametric in `FT` for parameters)
- Julia can infer types; use annotations primarily for **multiple dispatch**
- Tracer names used as type parameters should be wrapped in `Val` (e.g. `edible_detritus_name :: Val{EN}`)

## Numeric Types

- **Never hardcode Float64** in kernels/tendencies: no `0.0`, `1.0`
- Use `zero(FT)`, `one(FT)`, `convert(FT, 1//2)`, or rational literals, taking `FT` from the
  biogeochemistry's type parameter (e.g. `::NPD{FT}) where FT`) or from `eltype(grid)`
- Use `on_architecture` for data transfers — never manual `Array()` / `CuArray()` calls

## Memory Efficiency

- Favor inline computations over allocating temporary memory
- Design solutions that work within the existing framework (auxiliary fields, `update_biogeochemical_state!`)

## Staggered Grid & Indexing

- Velocities (including sinking velocities) live at cell faces, tracers at cell centers
- Take care of location when computing fluxes (e.g. sinking into sediments at the bottom face)
- **Always use 3D indexing** for fields (`field[i, j, k]`), including box models and columns
