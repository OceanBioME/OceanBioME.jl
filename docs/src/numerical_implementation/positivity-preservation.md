# [Positivity Preservation](@id pos-preservation)

It is common in BGC models to behave badly if any tracers go bellow zero, analytically this is fine because they (usually) can not get below zero, and it is unphysical for the concentration of something to be negative. Issues arise when the inaccuracy in numerical integration making some tracer become negative, usually leading to explosions (e.g. ``\exp(-C) \to \inf``), or bounds errors (e.g. ``\log(C)``). Essentially this occurs when the local error in the numerical scheme is sufficiently large that more than the available amount of tracer is consumed.

There exists a set of numerical schemes which overcome this and guarantee positivity (provided a positivity preserving advection scheme and well-behaved diffusion scheme are used) but are complex to implement. Although we may do this in the future we have not yet done so in the meantime have provided some utilities which maintain positivity. The simplest option is to reset any negative tracers to zero, but this causes the model to gain mass. A slightly more complex version is to increase negative tracers to zero and remove the difference from other tracers with in the same conserved system.

We have found this to be a satisfactory solution (when balanced against using much smaller time steps), as it tends to cause only a small and local transient change to the solution.

See model components page for usage.
## Conservative negative-tracer reservoirs

`ReservoirNegativeTracers` transfers negative concentrations into signed
reservoirs, one per eligible conservation group. For an elemental group with
tracer concentrations ``T_i`` and coefficients ``a_i``, the local correction
preserves ``\sum_i a_i T_i + R``, where ``R \leq 0`` is the reservoir. Deficits
can subsequently move through the domain and be repaid by a designated
nonnegative inorganic tracer. Including the reservoir in inventory diagnostics
is essential for measuring conservation.

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
set!(model.tracers.N, (x, y, z) -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
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
