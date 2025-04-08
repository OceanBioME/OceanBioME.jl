# [Sediment](@id sediment)

Sediment models can be added to biogeochemical models and evolve their own `biogeochemistry`.
They may have `fields` which can be prognostically integrated by the `timestepper`, and can store the value of tracers in the fluid domain as `tracked_fields` (for example sediments).
The `tracked_fields` may be `reqired_tracers`, or can enter the sediment model as `sinking_fluxes` where the `tracked_field` value will be:
```math
-w\frac{\partal C}{\partial z}.
```

```julia
struct BiogeochemicalSediment{BC, TS, CL, GR, SF, TF, BI} <: AbstractModel{TS}
    biogeochemistry :: BC
        timestepper :: TS
              clock :: CL
               grid :: GR
             fields :: SF
     tracked_fields :: TF
     bottom_indices :: BI
end
```

Please see the `InstantRemineralisationSediment` code as an example of how to implement a simple sediment biochemistry.

NB: `BiogeochemicalSediment` are currently only valid for flat bottom grids, or grids with flat cell bottoms (i.e. `PartialCellBottom` or `GridFittedBottom`), and will not be correct for future Oceananigans immersed boundaries such as cut-cell grids.