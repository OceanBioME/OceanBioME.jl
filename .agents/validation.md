# Validation Cases / Reproducing Published Results

Validation scripts live in `validation/` (e.g. `validation/LOBSTER/`, `validation/PISCES/`,
`validation/carbon_chemistry.jl`). They are not run on CI.

## 1. Parameter Extraction

- Use only the paper(s) supplied — never fill in parameters or references from memory; list
  anything missing and ask
- Extract: biogeochemical model structure and parameter values, domain and resolution, physical
  forcing (mixing / mixed-layer depth, surface PAR, temperature, salinity), initial conditions,
  boundary conditions (air–sea gas exchange, sediments, sinking), duration
- Check parameter tables, appendices and figure captions
- Record units for every parameter and convert to OceanBioME's (rates in 1/s, concentrations in
  mmol/m³ of the model currency); note sign conventions (Oceananigans `z` is negative
  downward, sinking velocities are negative, check air–sea flux sign)

## 2. Setup Verification (before long runs)

- Plot forcing functions (PAR, diffusivity, temperature) against time and depth
- Check grid extents and that `z` covers the intended depth
- After `set!`, check `minimum`/`maximum` of each tracer are physical and non-negative

## 3. Short Test Runs

- A few time steps on CPU at low resolution (or a `BoxModel` first)
- Choose `Δt` with `column_advection_timescale` / `sinking_advection_timescale`
- Check for NaNs and negative tracers; consider `ScaleNegativeTracers`
- Check conservation of the currency (e.g. total nitrogen) with closed boundaries
- Then run on GPU to catch GPU-specific issues

## 4. Progressive Validation

- Run a season or a year and visualize: bloom timing, nutrient drawdown, export, air–sea flux
- Compare with early-time figures in the paper

## 5. Comparison to Paper Figures

- Match figure format, units, axis ranges and time snapshots
- Compute the same diagnostics as the paper (e.g. integrated primary production, export flux)

## 6. Common Issues

- **NaNs / negative tracers**: `Δt` too large for fast rates or sinking, no positivity preservation
- **Nothing happening**: PAR zero (wrong `z` sign), tracers not set, biogeochemistry not passed to the model
- **Unit mismatches**: paper rates per day vs OceanBioME per second; mmol vs μmol; per kg vs per m³
- **GPU issues**: forcing functions must be `@inline`, type-stable and use `ifelse`
