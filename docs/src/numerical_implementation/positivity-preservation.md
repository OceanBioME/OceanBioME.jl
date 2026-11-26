# [Positivity Preservation](@id pos-preservation)

It is common in BGC models to behave badly if any tracers go bellow zero, analytically this is fine because they (usually) can not get below zero, and it is unphysical for the concentration of something to be negative. Issues arise when the inaccuracy in numerical integration making some tracer become negative, usually leading to explosions (e.g. ``\exp(-C) \to \inf``), or bounds errors (e.g. ``\log(C)``). Essentially this occurs when the local error in the numerical scheme is sufficiently large that more than the available amount of tracer is consumed.

There exists a set of numerical schemes which overcome this and guarantee positivity (provided a positivity preserving advection scheme and well-behaved diffusion scheme are used) but are complex to implement. Although we may do this in the future we have not yet done so in the meantime have provided some utilities which maintain positivity. 
OceanBioME.jl provides three modifiers for correcting negative tracer concentrations,
described below.

See [Utilities](@ref utils) for usage.

## ClipNegativeTracers

`ClipNegativeTracers` sets negative concentrations to zero without redistributing
the difference, increasing the affected tracer inventories.

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

## AntiNegativeTracers

`AntiNegativeTracers` sets negative tracer concentrations to zero while
conserving elemental inventories. Each selected conservation group (e.g.
nitrogen or carbon) receives one additional **anti-tracer**. Its value becomes
negative when the correction cannot be fully annihilated using a selected tracer.
Anti-tracers can be transported between grid cells, allowing annihilation
wherever the selected tracer is available.

For tracer concentrations ``T_i``, stoichiometric weights ``a_{gi}``, and the
anti-tracer value ``A_g \leq 0``, the correction preserves
``\sum_i a_{gi} T_i + A_g`` within each grid cell. Let
``T_i^+ = \max(T_i, 0)`` and let ``b_g`` denote the existing tracer selected
for annihilation in group ``g`` (specified by `annihilation_tracers`). Each
such tracer belongs exclusively to its own group. The correction applied
within each grid cell computes

```math
\begin{aligned}
\widetilde{A}_g &= A_g + \sum_i a_{gi} (T_i - T_i^+), \\
q_g &= T_{b_g}^+ + \widetilde{A}_g / a_{g b_g}, \\
T_{b_g}' &= \max(q_g, 0), \qquad A_g' = a_{g b_g} \min(q_g, 0).
\end{aligned}
```

Here, ``A_g`` is the anti-tracer value before correction, ``\widetilde{A}_g``
is its value after incorporating negative tracer concentrations but before
annihilation, and ``A_g'`` is its updated value after annihilation. Primes
denote updated values.

All other concentrations become ``T_i^+``.

The anti-tracers are ordinary Oceananigans tracers with zero biological
tendencies. Oceananigans handles their transport and boundary conditions.
The following example uses tracer-specific diffusion to illustrate anti-tracer
transport: diffusion of the original nutrient tracer is disabled, while
the anti-tracer diffuses between grid cells. Both tracers can undergo advection
and diffusion in general simulations. The simple abiotic configuration
isolates the correction from biological production and consumption.

```@example anti_tracer_diffusion
using OceanBioME, Oceananigans, CairoMakie
using Oceananigans.TimeSteppers: update_state!

grid = RectilinearGrid(size = 5, x = (0, 5), topology = (Bounded, Flat, Flat))

base = NPZD(grid; plankton = Abiotic(), detritus = InstantRemineralisationDetritus(),
            light_attenuation = nothing)
modifier = AntiNegativeTracers(base.underlying_biogeochemistry;
                              annihilation_tracers = (nitrogen = :N,))
bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)

closure = HorizontalScalarDiffusivity(ν = 0, κ = (N = 0, nitrogen_anti_tracer = 0.5))
model = NonhydrostaticModel(grid; biogeochemistry = bgc, closure,
                            momentum_advection = nothing,
                            tracer_advection = (N = Centered(), nitrogen_anti_tracer = nothing))

# Negative N initially occupies cell 1; cell 2 has no N; cells 3+ have N.
set!(model.tracers.N, x -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
set!(model.tracers.nitrogen_anti_tracer, 0)

profile(field) = vec(Array(interior(field)))
N_before = profile(model.tracers.N)
A_before = profile(model.tracers.nitrogen_anti_tracer)

update_state!(model) # N becomes nonnegative; the negative value is stored in Anti-N.
N_initial = profile(model.tracers.N)
A_initial = profile(model.tracers.nitrogen_anti_tracer)

simulation = Simulation(model; Δt = 0.1, stop_iteration = 60)
run!(simulation)
N_6 = profile(model.tracers.N)
A_6 = profile(model.tracers.nitrogen_anti_tracer)

simulation.stop_iteration = 500
run!(simulation)
N_final = profile(model.tracers.N)
A_final = profile(model.tracers.nitrogen_anti_tracer)
@assert isapprox(sum(N_before) + sum(A_before), sum(N_initial) + sum(A_initial); atol = 1e-9)
@assert isapprox(sum(N_initial) + sum(A_initial), sum(N_final) + sum(A_final); atol = 1e-9)
@assert sum(abs, A_final) < 1e-4 # Almost all Anti-N has been annihilated.
nothing
```

The saved concentrations are plotted before correction, immediately after
correction, and at two subsequent times.

```@example anti_tracer_diffusion
# Show N and Anti-N in every grid cell, before and after correction and transport.
cells = 1:5
stages = ((N_before, A_before), (N_initial, A_initial),
          (N_6, A_6), (N_final, A_final))
stage_titles = ("Before correction", "After correction (t = 0)",
                "Diffusion and annihilation (t = 6)", "Complete annihilation (t = 50)")

fig = Figure(size = (900, 950))
axes = Axis[]
color_map = cgrad([:darkorange, :white, :dodgerblue], [0, 1/3, 1])
for (stage, (N, A)) in enumerate(stages)
    ax = Axis(fig[stage, 1]; title = stage_titles[stage],
              yticks = ([1, 2], ["Anti-N", "N"]), xticks = cells,
              xlabel = stage == 4 ? "Grid cell" : "",
              xgridvisible = false, ygridvisible = false)
    push!(axes, ax)

    values = hcat(A[cells], N[cells]) # Anti-N below N on the same grid-cell axis.
    heatmap!(ax, cells, 1:2, values; colormap = color_map,
             colorrange = (-0.5, 1.0))
    for j in 1:2, i in eachindex(cells)
        value = values[i, j]
        label = abs(value) < 0.005 ? "0" : string(round(value; digits = 2))
        text!(ax, label; position = (cells[i], j), align = (:center, :center),
              color = abs(value) > 0.35 ? :white : :black, fontsize = 15)
    end
    vlines!(ax, 1.5:1:4.5; color = (:black, 0.2), linewidth = 1)
    hlines!(ax, [1.5]; color = (:black, 0.2), linewidth = 1)
    xlims!(ax, 0.5, 5.5)
    ylims!(ax, 0.5, 2.5)
    stage < 4 && hidexdecorations!(ax)
end
linkxaxes!(axes...)
Colorbar(fig[1:4, 2]; colormap = color_map, limits = (-0.5, 1.0),
         label = "Tracer value")
fig
```

The first stage shows the prescribed negative N concentration in cell 1, with
Anti-N initially zero. The correction sets N to zero there and stores its
negative value in Anti-N. Over time, Anti-N diffuses through cells without N;
where N is available, annihilation reduces its concentration. By ``t = 50``, all Anti-N has been annihilated. N itself does not diffuse in this example.
The figure shows all five grid cells. The combined inventory remains constant
under the no-flux boundaries. This example uses an explicit diffusion coefficient ``\kappa = 0.5`` and time step
``\Delta t = 0.1`` on a grid with ``\Delta x = 1``; larger diffusivities
require correspondingly smaller explicit time steps. In a three-dimensional
model, supply the desired tracer-specific advection and diffusivities for **all**
tracers when constructing the Oceananigans model.

### Conservation groups and annihilation

Each conservation group requires a tracer selected for annihilation that belongs
exclusively to that group. Tracers shared by multiple conservation groups are
corrected consistently.

When constructed from a biogeochemical model, `AntiNegativeTracers` selects
conservation groups with nonnegative stoichiometric weights. Groups with
negative weights (such as some oxygen groups) are excluded, while tracers
absent from the model's `conserved_tracers` groups retain their original
treatment. Conservation is guaranteed for the selected groups, including their
anti-tracer values.
