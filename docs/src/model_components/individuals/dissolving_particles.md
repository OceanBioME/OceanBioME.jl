# [Dissolving particles](@id dissolving-particles)

Dissolving particles are active particles which sink, dissolve, and release material into tracer fields.
They were written for modelling mineral-based ocean alkalinity enhancement (OAE), where particles of a mineral such as
lime (Ca(OH)₂), calcite (CaCO₃), or olivine (Mg₂SiO₄) are added to the ocean and increase the alkalinity as they dissolve, and to model the evolution of secondary mineral precipitates which can form in response to OAE.

The only property tracked by each particle is its radius, ``r`` (m), which sets how fast the particle dissolves and how fast it sinks.
An example of particles sinking and dissolving in two-dimensional turbulence is given in the [examples](@ref dissolving-particles-example).

You can set up the particle biogeochemistry with [`DissolvingMineral`](@ref):
```jldoctest
using OceanBioME

lime = DissolvingMineral()

# output
DissolvingMineral{Float64} biogeochemistry tracking particle `r`adius
├── density: 2211.0 kg/m³
├── molar mass: 0.074093 kg/mol
├── released tracers: (Alk = 2.0,) mol/mol
├── dissolution rate: 1.0e-8 m/s
└── sinking speed: StokesSinking{Float64}(2211.0, 1025.0, 0.00108, 9.81) m/s
```
which can be put into `BiogeochemicalParticles`, or you can directly create the particles with [`DissolvingParticles`](@ref):
```jldoctest
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1));
particles = DissolvingParticles(10; grid)

# output
10 BiogeochemicalParticles with DissolvingMineral{Float64} biogeochemistry:
├── fields: (:r,)
└── coupled tracers: (:Alk,)
```

The particles can then be added to any biogeochemical model which has the tracers they release,
for example a carbonate system without any biology:
```julia
biogeochemistry = NutrientsPlanktonDetritus(grid; inorganic_carbon = CarbonateSystem(), particles)
```

## Model equations

The radius of each particle shrinks at the dissolution rate, ``k(r)`` (m/s):

```math
\frac{dr}{dt} = -k(r),
```

so the volume of each particle, ``V = 4\pi r^3/3``, decreases at ``4\pi r^2 k(r)``.
When the dissolution is controlled by reactions on the particle surface ``k`` is a constant
(the default), but it can be any function of radius. For example, when dissolution is limited by diffusion
away from the particle surface ``k \propto 1/r``.

As the particles dissolve they release each tracer at:

```math
\frac{dX}{dt} = \frac{4\pi r^2 k(r) \nu_X}{V_m},
```

where ``\nu_X`` is the number of moles of tracer ``X`` released per mole of mineral (`released_tracers`), and
``V_m`` is the molar volume of the mineral (`molar_mass / density`).
For example, calcium hydroxide releases two moles of alkalinity per mole (`released_tracers = (; Alk = 2)`),
and calcium carbonate releases two moles of alkalinity and one mole of dissolved inorganic carbon
(`released_tracers = (; Alk = 2, DIC = 1)`).

The released tracers go into the grid cell nearest to the particle, and each particle can represent many real
particles by setting its `scalefactor`.

### Sinking

When the particles are moved by `LagrangianAdvection` (the default) they move with the water, and sink relative to it at `sinking_speed`,
which can be a number or a function of radius. By default they sink at the speed given by Stokes' law
(see [`StokesSinking`](@ref)):

```math
w = \frac{2(\rho_p - \rho)gr^2}{9\mu},
```

where ``\rho_p`` is the density of the particle, ``\rho`` the density of the water, ``g`` the gravitational acceleration,
and ``\mu`` the dynamic viscosity of the water. This is only valid for small particles (with a Reynolds number,
``2rw\rho/\mu``, less than about one), so a different `sinking_speed` should be given for larger particles.

Particles which reach the bottom of the domain bounce off it, so they stay within one timestep's sinking of the bottom
and continue to dissolve into the bottom grid cells.

### Parameter variable names

| Symbol      | Variable name      | Units                  |
|-------------|--------------------|------------------------|
| ``\rho_p``  | `density`          | kg / m³                |
| -           | `molar_mass`       | kg / mol               |
| ``\nu_X``   | `released_tracers` | mol / mol              |
| ``k``       | `dissolution_rate` | m / s                  |
| ``w``       | `sinking_speed`    | m / s                  |

The default `density`, `molar_mass`, and `released_tracers` are for calcium hydroxide. The default dissolution rate is only illustrative
(a 100 μm particle dissolves in about 3 hours) and should be set for the mineral and conditions being modelled.

## Model conservations

The total of each released tracer, in the water and in the particles, is conserved. The amount of tracer ``X`` in each particle is:

```math
X_\text{particle} = \frac{4\pi r^3 \nu_X}{3V_m}.
```

Note that because the particle radius and the tracers are stepped forward differently, this is conserved up to a first order timestepping error. When the particles fully dissolve slightly too much tracer is released, about ``0.5k\Delta t / r_0`` of the total with the default `RungeKutta3` timestepper (or ``3k\Delta t / r_0`` with `QuasiAdamsBashforth2`). As a consequence, the timestep should be much shorter than the time it takes the particles to dissolve, ``r_0/k`` to maintain material conservation.
