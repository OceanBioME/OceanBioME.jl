# MARBL Fortran comparison

These compare the OceanBioME MARBL components against the **standalone MARBL Fortran driver's**
baseline output, at two levels:

- **`compare_marbl.jl`** — individual rate terms and reconstructed tendencies, section by section
  (per-rate terms and the grazing network, the real tendency methods, assembled `J_<tracer>`,
  DOM/POM/PFe, carbonate and calcite, O₂/N redox, +cocco, the iron cycle, sediments).
- **`test_stepped_marbl.jl`** — the strongest check: each MARBL column is built on its own
  bathymetry-deep grid and stepped through the real Oceananigans timestepper, and `Gⁿ` is compared
  against MARBL's `J_<tracer>` at every interior cell, sediment coupling included, for four
  configurations (base CESM2.1, +cocco, and two general N-plankton sets).

## Why these are not in `runtests.jl`

They read MARBL's baseline `.nc` files (~8 MB, under `marbl-info3/MARBL/tests/input_files/` and
`marbl-info3/draft/`), which are not part of the package, and they need `NCDatasets`. They are
developer gates, run deliberately.

## Running them

```
julia --check-bounds=yes --project=test/marbl test/marbl/runtests.jl
```

That runs both, each in its own process, and prints one summary at the end. Either can be run alone:

```
julia --check-bounds=yes --project=test/marbl test/marbl/compare_marbl.jl
julia --check-bounds=yes --project=test/marbl test/marbl/test_stepped_marbl.jl
```

**Read the summary, not the running output.** Each takes tens of minutes and prints as it goes, so a
log inspected mid-run will show you the previous section's numbers — or, if you have just edited
`src`, the *previous run's* output while this one recompiles. `runtests.jl` reports only on
completion for exactly this reason.

If the baselines are missing both skip with a warning rather than failing, so the directory is safe
to run in a checkout that does not have them.

`--check-bounds=yes` matters: several of these checks read fields near the sea floor and the domain
edges, where an out-of-bounds read would otherwise pass silently.

## `marbl_names.jl`

The comparison logic was written against the standalone prototype, whose types carried MARBL-specific
names. The components are now general OceanBioME ones with MARBL as a parameter set, so
`marbl_names.jl` is the single place the two vocabularies meet — `ManyPhytoZoo as MARBLPlankton` and
so on. It also carries the handful of signatures that changed in the port, and two things worth
knowing when reading a failure:

- **Temperature is not a biogeochemical tracer.** The models here declare `T` themselves
  (`PHYSICS_TRACERS`) and set it from the baseline. If a harness forgets, every temperature-dependent
  rate silently evaluates at 0 °C and *every* autotroph tendency fails by order unity.
- **Tracer names carry subscripts** (`spCaCO₃`), MARBL's baseline variables are ASCII (`spCaCO3`).
  `ascii_name` converts when reading a baseline. Passing an ASCII `Val` *into* the model is worse
  than a lookup error: no method matches, so the abstraction's zero default is returned and the
  comparison reports `ours = 0`.

## Layout

| file | |
|---|---|
| `runtests.jl` | runs both and reports once |
| `marbl_baselines.jl` | where the baselines are, how to read them, the MARBL-matched carbon chemistry |
| `marbl_names.jl` | the name aliases and the signatures that changed in the port |
| `general_configs.jl` | the general N-plankton configurations |
| `compare_marbl.jl` | rate terms, section by section |
| `test_stepped_marbl.jl` | assembled tendencies through the timestepper |

These began as copies of the equivalents in `marbl-info3/draft/`, which is the prototype the
components were ported from, and are now maintained here.
