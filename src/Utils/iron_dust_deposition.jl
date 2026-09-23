using Oceananigans.BoundaryConditions: getbc
using Oceananigans.Forcings: Forcing
using Oceananigans.Grids: znode, Center

import Adapt: adapt_structure, adapt
import Base: summary, show

"""
    IronDustDeposition(dust_deposition;
                       iron_mass_fraction = 0.035,
                       molar_mass = 55.845e-6,
                       dissolution_length = 400.0,
                       hard_fraction = 0.98,
                       hard_dissolution_length = 1.2e6)

Implicit model of iron dissolution from sinking dust particles. Dust deposited at the surface
sinks through the water column and dissolves, releasing iron with an exponential decay profile:

    source(z) = F₀ × fᵢ / M × ((1 - γ) / λₛ × exp(z/λₛ)  +  γ / λₕ × exp(z/λₕ))

where `F₀` is the surface dust flux (kg/m²/s), `fᵢ` is the iron mass fraction, `M` is the
molar mass, `γ` is the hard fraction, and `λₛ`, `λₕ` are the dissolution lengths.

Defaults match MARBL's two-component formulation where 98% of dust is mineral-ballasted and
barely dissolves. Set `hard_fraction = 0` for a single-exponential model.

`dust_deposition` can be a number, a `Field`, a `FieldTimeSeries`, or anything `getbc` accepts.

Use `IronDustDepositionForcing` to get a ready-made `Forcing` for the `Fe` tracer.

Keyword Arguments
=================

- `iron_mass_fraction`: mass fraction of iron in dust (default 0.035, i.e. 3.5%)
- `molar_mass`: molar mass of iron in kg/mmol (default 55.845e-6)
- `dissolution_length`: e-folding depth for the (soft) fraction in meters (default 400.0)
- `hard_fraction`: fraction of dust entering a slowly dissolving subclass (default 0.98)
- `hard_dissolution_length`: e-folding depth for the hard fraction in meters (default 1.2e6)
"""
struct IronDustDeposition{DF, FT} <: Function
             dust_deposition :: DF
          iron_mass_fraction :: FT
                  molar_mass :: FT
        dissolution_length :: FT
               hard_fraction :: FT
    hard_dissolution_length :: FT
end

IronDustDeposition(dust_deposition;
                   iron_mass_fraction = 0.035,
                   molar_mass = 55.845e-6,
                   dissolution_length = 400.0,
                   hard_fraction = 0.98,
                   hard_dissolution_length = 1.2e6) =
    IronDustDeposition(dust_deposition,
                       iron_mass_fraction,
                       molar_mass,
                       dissolution_length,
                       hard_fraction,
                       hard_dissolution_length)

adapt_structure(to, d::IronDustDeposition) =
    IronDustDeposition(adapt(to, d.dust_deposition),
                       adapt(to, d.iron_mass_fraction),
                       adapt(to, d.molar_mass),
                       adapt(to, d.dissolution_length),
                       adapt(to, d.hard_fraction),
                       adapt(to, d.hard_dissolution_length))

summary(::IronDustDeposition) = "IronDustDeposition"
function show(io::IO, d::IronDustDeposition)
    print(io, summary(d), "\n",
          "├── Dust deposition: ", summary(d.dust_deposition), "\n",
          "├── Iron mass fraction: ", d.iron_mass_fraction, "\n",
          "├── Molar mass: ", d.molar_mass, " kg/mmol\n",
          "├── Dissolution length: ", d.dissolution_length, " m")
    if d.hard_fraction > 0
        print(io, "\n├── Hard fraction: ", d.hard_fraction,
                  "\n└── Hard dissolution length: ", d.hard_dissolution_length, " m")
    end
end

@inline function (d::IronDustDeposition)(i, j, k, grid, clock, model_fields)
    F₀ = getbc(d.dust_deposition, i, j, grid, clock, model_fields)
    z  = znode(i, j, k, grid, Center(), Center(), Center())
    γ  = d.hard_fraction
    λₛ = d.dissolution_length
    dust_to_iron = F₀ * d.iron_mass_fraction / d.molar_mass
    soft = (1 - γ) / λₛ * exp(z / λₛ)
    hard = γ / d.hard_dissolution_length * exp(z / d.hard_dissolution_length)
    return dust_to_iron * (soft + hard)
end

"""
    IronDustDepositionForcing(dust_deposition; kwargs...)

Returns a discrete `Forcing` for the `Fe` tracer that models iron dissolution from sinking dust
using MARBL's two-component formulation by default.

`dust_deposition` is the surface dust deposition rate in kg/m²/s. It can be a constant, a 2D
`Field`, a `FieldTimeSeries`, or a function.

```julia
forcing = (; Fe = IronDustDepositionForcing(dust_field))
model = NonhydrostaticModel(grid; biogeochemistry, forcing)
```
"""
IronDustDepositionForcing(dust_deposition; kwargs...) =
    Forcing(IronDustDeposition(dust_deposition; kwargs...); discrete_form = true)
