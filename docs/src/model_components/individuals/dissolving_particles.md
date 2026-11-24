# [Dissolving particles](@id dissolving-particles)

Dissolving particles are active particles which sink, dissolve, and release material into tracer fields.
They were written for mineral-based ocean alkalinity enhancement (OAE), where particles of a mineral such as
lime (Ca(OH)₂), calcite (CaCO₃), or olivine (Mg₂SiO₄) increase the alkalinity as they dissolve, and for
the secondary mineral precipitates which can form in response to OAE.

Each particle tracks its radius, ``r`` (m), which sets how fast it dissolves and sinks.
An example of particles sinking and dissolving in two-dimensional turbulence is given in the
[examples](@ref dissolving-particles-example).

The particle biogeochemistry is set up with [`DissolvingMineral`](@ref):
```@example dissolving_particles
using OceanBioME, Oceananigans

mineral = DissolvingMineral()
```
which can be put into `BiogeochemicalParticles`, or the particles can be created directly with
[`DissolvingParticles`](@ref):
```@example dissolving_particles
grid = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1))

particles = DissolvingParticles(10; grid)
```

The particles can then be added to any biogeochemical model which has the tracers they release,
for example a carbonate system without any biology:
```@example dissolving_particles
biogeochemistry = NutrientsPlanktonDetritus(grid; inorganic_carbon = CarbonateSystem(), particles)
```

## Model equations

The radius of each particle shrinks at the dissolution rate, ``k`` (m/s):

```math
\frac{dr}{dt} = -k,
```

so the volume of each particle, ``V = 4\pi r^3/3``, decreases at ``4\pi r^2 k``.
When the dissolution is controlled by reactions on the particle surface ``k`` is a constant
(the default), and when it is limited by diffusion away from the particle surface ``k \propto 1/r``.

As the particles dissolve they release each tracer at:

```math
\frac{dX}{dt} = \frac{4\pi r^2 k \theta_X}{V_m},
```

where ``\theta_X`` is the number of moles of tracer ``X`` released per mole of mineral (`released_tracers`), and
``V_m`` is the molar volume of the mineral (`molar_mass / density`).
For example, calcium hydroxide releases two moles of alkalinity per mole (`released_tracers = (; Alk = 2)`),
and calcium carbonate releases two moles of alkalinity and one mole of dissolved inorganic carbon
(`released_tracers = (; Alk = 2, DIC = 1)`).

The released tracers go into the grid cell nearest to the particle, and each particle can represent many real
particles by setting its `scalefactors`.

### Sinking

When the particles are moved by `LagrangianAdvection` (the default) they move with the water, and also
move vertically relative to it at `sinking_speed`, which is positive upwards, so it is negative for sinking
particles. By default it is given by Stokes' law ([`StokesSinking`](@ref)):

```math
w = \frac{2(\rho - \rho_p)gr^2}{9\mu},
```

where ``\rho_p`` is the density of the particle, ``\rho`` the density of the water, ``g`` the gravitational acceleration,
and ``\mu`` the dynamic viscosity of the water. This is only valid for small particles, with a Reynolds number
``2r|w|\rho/\mu`` less than about one.

The `sinking_speed` of each particle is stored in the particle field of the same name, which is updated every
time step. Particles which reach the bottom of the domain bounce off it, so they stay within one time step's
sinking of the bottom and continue to dissolve into the bottom grid cells.

### Dependence on particle and model fields

`dissolution_rate` and `sinking_speed` can be numbers, or functions called as
`f(mineral, particle_fields, tracer_fields)`, where `particle_fields` is a `NamedTuple` of the particle's
`r` and `sinking_speed`, and `tracer_fields` is a `NamedTuple` of the model fields listed in `required_fields`
(tracers, velocities, or auxiliary fields) at the particle position.
For example, a dissolution rate which increases with temperature:
```@example dissolving_particles
warm_dissolution(mineral, particle_fields, tracer_fields) = 1e-8 * (1 + 0.05 * (tracer_fields.T - 20))

mineral = DissolvingMineral(required_fields = (:T, ), dissolution_rate = warm_dissolution)
```

The `water_density` and `dynamic_viscosity` of [`StokesSinking`](@ref) can be functions called in the same way,
for example for particles which sink faster as the water warms and becomes less viscous:
```@example dissolving_particles
viscosity(mineral, particle_fields, tracer_fields) = 1.8e-3 / (1 + 0.03 * tracer_fields.T)

stokes = StokesSinking(particle_density = 2211, dynamic_viscosity = viscosity)

mineral = DissolvingMineral(required_fields = (:T, ), sinking_speed = stokes)
```

### Parameter variable names

| Symbol      | Variable name      | Units                  |
|-------------|--------------------|------------------------|
| ``\rho_p``  | `density`          | kg / m³                |
| -           | `molar_mass`       | kg / mol               |
| ``\theta_X`` | `released_tracers` | mol / mol              |
| ``k``       | `dissolution_rate` | m / s                  |
| ``w``       | `sinking_speed`    | m / s                  |
| -           | `required_fields`  | -                      |

The default `density`, `molar_mass`, and `released_tracers` are for calcium hydroxide.
The default dissolution rate is only illustrative (a 100 μm particle dissolves in about 3 hours).

## Model conservations

The total of each released tracer, in the water and in the particles, is conserved. The amount of tracer ``X`` in each particle is:

```math
X_\text{particle} = \frac{4\pi r^3 \theta_X}{3V_m}.
```

Because the particle radius and the tracers are stepped forward differently, this is conserved up to a first
order time stepping error. When the particles fully dissolve slightly too much tracer is released, about
``0.5k\Delta t / r_0`` of the total with the default `RungeKutta3` time stepper (or ``3k\Delta t / r_0`` with
`QuasiAdamsBashforth2`), so the time step should be much shorter than the time it takes the particles to
dissolve, ``r_0/k``.
