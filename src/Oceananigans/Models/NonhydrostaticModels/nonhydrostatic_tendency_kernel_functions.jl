"""
Return the tendency for a auxiliary field with index `tracer_index`
at grid point `i, j, k`.

The tendency is called ``G_a`` and defined via

```math
dc_t a = G_a ,
```

where `a = A[tracer_index]`.

The arguments `velocities` and `tracers` are `NamedTuple`s with the three
velocity components and tracer fields where applicable.
`forcings` is a named tuple of forcing functions.

`clock` keeps track of `clock.time` and `clock.iteration`.
"""
@inline function auxiliary_tendency(i, j, k, grid,
                                 velocities,
                                 tracers,
                                 auxiliary_fields,
                                 forcing,
                                 clock)

    model_fields = merge(velocities, tracers, auxiliary_fields)

    return forcing(i, j, k, grid, clock, model_fields)
end