using Oceananigans.BoundaryConditions: FluxBoundaryCondition

import Adapt: adapt_structure, adapt
import Base: summary, show

"""
    IronDustDeposition(dust_iron_flux; solubility = 0.01)

A callable struct for use as a discrete-form `FluxBoundaryCondition` on the `Fe` tracer, following
the MITgcm DIC formulation for aeolian iron deposition.

The iron flux into the surface cell is `solubility × dust_iron_flux`, where `dust_iron_flux` is the
total iron flux in dust (in the same concentration units as `Fe`, e.g. mmol Fe/m²/s) and `solubility`
is the fraction of deposited iron that dissolves (default 1%, following Jickells et al., 2005).

`dust_iron_flux` can be a number, a function of `(i, j, grid, clock, model_fields)`, a `Field`, or
a `FieldTimeSeries`.

Returns a `FluxBoundaryCondition` when called with
[`IronDustDepositionBoundaryCondition`](@ref).

Keyword Arguments
=================

- `solubility`: fraction of deposited iron that dissolves (default 0.01)
"""
struct IronDustDeposition{DF, FT} <: Function
    dust_iron_flux :: DF
        solubility :: FT
end

IronDustDeposition(dust_iron_flux; solubility = 0.01) =
    IronDustDeposition(dust_iron_flux, solubility)

adapt_structure(to, d::IronDustDeposition) =
    IronDustDeposition(adapt(to, d.dust_iron_flux), adapt(to, d.solubility))

summary(::IronDustDeposition) = "IronDustDeposition"
show(io::IO, d::IronDustDeposition) =
    print(io, summary(d), "\n",
          "├── Dust iron flux: ", summary(d.dust_iron_flux), "\n",
          "└── Solubility: ", d.solubility)

@inline (d::IronDustDeposition)(i, j, grid, clock, model_fields) =
    - d.solubility * surface_value(d.dust_iron_flux, i, j, grid, clock, model_fields)

@inline surface_value(f::Number, i, j, grid, clock, model_fields) = f
@inline surface_value(f::Function, i, j, grid, clock, model_fields) = f(i, j, grid, clock, model_fields)
@inline surface_value(f::AbstractArray{<:Any, 2}, i, j, grid, clock, model_fields) = @inbounds f[i, j]

"""
    IronDustDepositionBoundaryCondition(dust_iron_flux; solubility = 0.01)

Returns a `FluxBoundaryCondition` for aeolian iron deposition on the `Fe` tracer.

`dust_iron_flux` is the total iron content of deposited dust (in the model's iron units per m²
per second, e.g. mmol Fe/m²/s). Only the fraction `solubility` dissolves and enters the ocean.

Usage:

```julia
iron_bc = IronDustDepositionBoundaryCondition(dust_iron_flux; solubility = 0.01)

# then pass to the model as:
# Fe = FieldBoundaryConditions(top = iron_bc)
```
"""
IronDustDepositionBoundaryCondition(dust_iron_flux; solubility = 0.01) =
    FluxBoundaryCondition(IronDustDeposition(dust_iron_flux; solubility); discrete_form = true)
