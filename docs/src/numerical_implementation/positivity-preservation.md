# [Positivity Preservation](@id pos-preservation)

It is common in BGC models to behave badly if any tracers go bellow zero, analytically this is fine because they (usually) can not get below zero, and it is unphysical for the concentration of something to be negative. Issues arise when the inaccuracy in numerical integration making some tracer become negative, usually leading to explosions (e.g. ``\exp(-C) \to \inf``), or bounds errors (e.g. ``\log(C)``). Essentially this occurs when the local error in the numerical scheme is sufficiently large that more than the available amount of tracer is consumed.

Positivity-preserving numerical schemes exist, but can be complex to implement.
OceanBioME.jl provides three modifiers for correcting negative tracer concentrations.
The simplest, [`ClipNegativeTracers`](@ref), resets negative tracers to zero,
but increases the affected tracer inventories. [`ScaleNegativeTracers`](@ref)
instead increases negative tracers to zero and removes the difference from other
tracers in the same conservation group, preserving its weighted total when
sufficient nonnegative inventory is available.

These corrections can introduce local, transient changes to the solution;
smaller time steps can reduce how frequently they are needed.

[`ReservoirNegativeTracers`](@ref) instead retains unresolved deficits in
additional prognostic reservoir tracers, allowing deficits to be transported
and subsequently repaid while conserving the selected elemental inventories.

See [Utilities](@ref utils) for usage.

## ClipNegativeTracers

`ClipNegativeTracers` sets negative concentrations to zero without redistributing
the difference.

## ScaleNegativeTracers

For a single conservation group, `ScaleNegativeTracers` sets
negative concentrations to zero and multiplies positive concentrations by
``t/p``, where ``t`` is the original weighted total and ``p`` is the weighted
total of positive concentrations. A negative group total cannot be preserved
using only nonnegative concentrations.

Tracers may belong to multiple conservation groups, such as nitrogen and carbon.
These groups must be passed to the same `ScaleNegativeTracers` modifier: separate
modifiers act sequentially and may change totals preserved by earlier corrections.
Each group needs at least one tracer that is not shared with another group.

## ReservoirNegativeTracers

`ReservoirNegativeTracers` transfers negative concentrations into signed
reservoirs, one per eligible conservation group. For an elemental group with
tracer concentrations ``T_i`` and coefficients ``a_i``, the local correction
preserves ``\sum_i a_i T_i + R``, where ``R \leq 0`` is the reservoir. Deficits
can subsequently move through the domain and be repaid by a designated
nonnegative inorganic tracer. Including the reservoir in inventory diagnostics
is essential for measuring conservation.

Let ``T_i^+ = \max(T_i, 0)`` and let ``b_g`` denote the repayment tracer for
group ``g``. Each repayment tracer belongs exclusively to its own group. The
correction computes

```math
\begin{aligned}
d_g &= R_g + \sum_i a_{gi} (T_i - T_i^+), \\
q_g &= T_{b_g}^+ + d_g / a_{g b_g}, \\
T_{b_g}' &= \max(q_g, 0), \qquad R_g' = a_{g b_g} \min(q_g, 0).
\end{aligned}
```

All other concentrations become ``T_i^+``. This preserves each selected
group's elemental inventory during the local correction.

The reservoirs are **ordinary Oceananigans prognostic tracers** with zero
biological tendencies. Oceananigans handles their transport and boundary
conditions. The following example uses its tracer-specific advection and
horizontal diffusion options: only the reservoir diffuses, and neither field
has a prescribed velocity. The simple abiotic configuration isolates the
correction from biological production and consumption.

```@example reservoir_diffusion
using OceanBioME, Oceananigans, CairoMakie
using Oceananigans.TimeSteppers: update_state!

grid = RectilinearGrid(size = 16, x = (0, 16), topology = (Bounded, Flat, Flat))

base = NPZD(grid; plankton = Abiotic(), detritus = InstantRemineralisationDetritus(),
            light_attenuation = nothing)
modifier = ReservoirNegativeTracers(base.underlying_biogeochemistry;
                                   repayment_tracers = (nitrogen = :N,))
bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)

closure = HorizontalScalarDiffusivity(ν = 0, κ = (N = 0, nitrogen_reservoir = 0.5))
model = NonhydrostaticModel(grid; biogeochemistry = bgc, closure,
                            momentum_advection = nothing,
                            tracer_advection = (N = Centered(), nitrogen_reservoir = nothing))

# Negative N initially occupies cell 1; cell 2 has no N; cells 3+ have N.
set!(model.tracers.N, x -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
update_state!(model) # N is nonnegative; the initial deficit is now in the reservoir.

profile(field) = vec(Array(interior(field)))
N_initial = profile(model.tracers.N)
R_initial = profile(model.tracers.nitrogen_reservoir)

simulation = Simulation(model; Δt = 0.1, stop_iteration = 12)
run!(simulation)

N_final = profile(model.tracers.N)
R_final = profile(model.tracers.nitrogen_reservoir)
@assert isapprox(sum(N_initial) + sum(R_initial), sum(N_final) + sum(R_final); atol = 1e-9)

x = 0.5:1:15.5
fig = Figure(size = (800, 390))
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "Concentration")
lines!(ax, x, N_initial; label = "N: initially", linestyle = :dash, color = :dodgerblue)
lines!(ax, x, N_final; label = "N: after diffusion", color = :dodgerblue)
lines!(ax, x, R_initial; label = "Reservoir: initially", linestyle = :dash, color = :darkorange)
lines!(ax, x, R_final; label = "Reservoir: after diffusion", color = :darkorange)
axislegend(ax; position = :rb)
fig
```

The deficit first enters a nutrient-free neighbouring cell as negative reservoir
inventory. Upon reaching cells with available N, it reduces N instead. The
combined inventory remains constant under the no-flux boundaries. This example
uses an explicit diffusion coefficient ``\kappa = 0.5`` and time step
``\Delta t = 0.1`` on a grid with ``\Delta x = 1``; larger diffusivities
require correspondingly smaller explicit time steps. In a three-dimensional
model, supply the desired tracer-specific advection and diffusivities for **all**
tracers when constructing the Oceananigans model.

### Conservation groups and repayment tracers

Each conservation group requires a repayment tracer that belongs exclusively to
that group. Tracers shared by multiple conservation groups are corrected
consistently.

When constructed from a biogeochemical model, `ReservoirNegativeTracers` selects
conservation groups with nonnegative stoichiometric coefficients. Groups with
signed coefficients (such as oxygen) are excluded, while tracers absent from the
model's `conserved_tracers` groups retain their original treatment. Conservation
is guaranteed for the selected groups, including their reservoirs.
