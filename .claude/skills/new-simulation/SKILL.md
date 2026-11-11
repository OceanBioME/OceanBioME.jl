---
name: new-simulation
description: Set up, run, and visualize a new OceanBioME simulation (box, column, or 2D/3D), with or without a reference paper
---

# New Simulation

## Step 1: Understand the Case

**If reproducing a paper:**
- Read the paper the user supplies and extract ALL parameters: biogeochemical model and parameter
  values, domain, resolution, physical forcing (mixing, PAR, temperature, salinity), initial
  conditions, boundary conditions (air–sea gas exchange, sediment), run duration
- Record units and sign conventions for each (e.g. `z` is negative downward in Oceananigans;
  sinking velocities are negative; air–sea fluxes: check the sign convention of the boundary condition)
- Don't fill gaps from memory — list missing parameters and ask

**If designing a new case:**
- Ask for the science goal, then clarify: which biogeochemical model (`NPZD`, `LOBSTER`,
  `PISCES`, a custom NPD configuration), light model, carbon chemistry / gas exchange, sediments,
  particles, dimensionality, and duration
- Identify what to diagnose (tracer profiles, fluxes, budgets via `Budget`, export)

## Step 2: Pick the Model Type

- **0D**: `BoxModel(; biogeochemistry, clock)` on a `BoxModelGrid()` with
  `PrescribedPhotosyntheticallyActiveRadiation` — see `examples/box.jl`
- **1D column**: `RectilinearGrid(topology = (Flat, Flat, Bounded), ...)` with
  `NonhydrostaticModel` or `HydrostaticFreeSurfaceModel`, prescribed diffusivity — see
  `examples/column.jl`
- **2D/3D with dynamics**: see `examples/eady.jl`; GPU recommended

## Step 3: Set Up Grid, Forcing and Initial Conditions

- Build the grid and check extents; for columns, check that `z` spans the intended depth
- Plot forcing functions (surface PAR, mixed-layer depth, diffusivity) against time before using them
- After `set!`, check `minimum`/`maximum` of each tracer make physical sense (concentrations in the
  model's units, e.g. mmol N/m³) and are non-negative
- Consider `ScaleNegativeTracers` in `modifiers` for long runs (see
  `docs/src/numerical_implementation/positivity-preservation.md`)

## Step 4: Short Test Run

- A few time steps on CPU at low resolution; choose `Δt` with `column_advection_timescale` /
  `sinking_advection_timescale` for sinking tracers
- Check for NaNs and negative tracers
- Check conservation of the currency (e.g. total nitrogen) over the run when there are no
  sources/sinks at the boundaries
- Check output files contain meaningful data

## Step 5: Progressive Validation

- Short simulation, visualize, check bloom timing, nutrient drawdown, export, air–sea flux magnitude
- If reproducing a paper, compare to early-time figures

## Step 6: Production Run and Comparison

- Full resolution / duration; diagnostics matching the science goal
- If reproducing a paper, match figure format, axis ranges, units, and time snapshots

## Visualization

Existing examples plot with `interior(...)` and explicit `times`/`z` arrays, which is fine. With
the Oceananigans Makie extension (loaded by `using CairoMakie` with Oceananigans) you can also plot
`Field`s directly — check the Oceananigans docs for the installed version before relying on it.

```julia
using OceanBioME, Oceananigans, CairoMakie
using Oceananigans.Units

P = FieldTimeSeries("column.jld2", "P")
times = P.times
x, y, z = nodes(P)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Time (days)", ylabel = "z (m)")
heatmap!(ax, times / days, z, interior(P, 1, 1, :, :)'; colormap = :batlow)
```

For box models, `FieldTimeSeries("box.jld2", "P")[1, 1, 1, :]` gives the time series
(see `examples/box.jl`).

Animations: use an `Observable` index and `@lift`, then `record(fig, "file.mp4", 1:length(times))`.

## Common Issues

- **NaNs / negative tracers**: `Δt` too large (especially with fast rates or sinking), missing
  positivity preservation, unstable ICs
- **Nothing happening**: PAR zero or wrong sign of `z`, tracers not initialised, biogeochemistry
  not passed to the model
- **Sinking tracers pile up or vanish at the bottom**: check `open_bottom` and sediment coupling
- **GPU issues**: forcing functions must be `@inline`, type-stable, and use `ifelse`

## Output

- Validation scripts: `validation/<case_name>/`
- New documented examples: `examples/` (see `.claude/rules/examples-rules.md` for registering them)
