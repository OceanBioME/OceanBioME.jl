---
paths:
  - src/**/*.jl
---

# Docstring Rules

## Format

OceanBioME does not depend on DocStringExtensions.jl. Write the signature(s) by hand:

~~~~
"""
    PhytoZoo([FT = Float64;] kwargs...)
    PhytoZoo(grid; phytoplankton_sinking_speed = 0, kwargs...)

One-paragraph description of what the component does and which tracers it adds.

Keyword Arguments
=================

- `maximum_grazing_rate`: maximum zooplankton grazing rate (1/s)
- `grazing_half_saturation`: grazing half-saturation constant (mmol N/m³)
"""
~~~~

- Give units for every dimensional keyword argument
- Cross-reference with `[`Name`](@ref)`
- Only cite papers the user has provided or that are already in `docs/oceanbiome.bib`;
  never add references from memory

## Prefer `jldoctest` over plain `julia` blocks

Plain code blocks (`` ```julia ``) are not tested and go stale. Doctests (`` ```jldoctest ``)
run in the CPU test suite (`doctest(OceanBioME)` at the end of `test/runtests.jl`).
Note `docs/make.jl` currently has `doctest = false`, so the docs build does **not** check them.

~~~~
"""
    CarbonChemistry(FT = Float64; kwargs...)

Example
=======

```jldoctest
julia> using OceanBioME

julia> carbon_chemistry = CarbonChemistry()
`CarbonChemistry` model which solves for pCO₂ and pH

julia> pCO₂ = carbon_chemistry(; DIC = 2000.0, Alk = 2000.0, T = 10.0, S = 35.0)
1308.1006995915372
```
"""
~~~~

(from `src/Models/CarbonChemistry/carbon_chemistry.jl`). Copy the exact output from running the
code — do not guess what `show` prints.

## Doctest Best Practices

- Always include expected output (either REPL-style `julia>` blocks or a `# output` section)
- Doctests should exercise `Base.show` to verify objects display correctly
- Use `filter = r"..."` for output that varies (paths, timings), as in `simple_multi_G.jl`
- Keep doctests minimal and CPU-only
- **Do NOT use boolean comparisons as the final line** (e.g. avoid `x ≈ 1.0` or `obj isa Type`)
