import Oceananigans.BoundaryConditions: getbc

"""
    PARFromShortwave(surface_shortwave; photosynthetic_fraction_of_shortwave = 0.43)

A surface photosynthetically active radiation (`PAR`) which is diagnosed as a fixed fraction of the
surface downwelling shortwave radiation, i.e.

```math
PAR_0 = f_{PAR} Q_{sw},
```

where ``f_{PAR}`` is `photosynthetic_fraction_of_shortwave`.

`surface_shortwave` may be anything the boundary condition `getbc` interface accepts (a number, a
`Field`, a `FieldTimeSeries`, ...), and is passed as the `surface_PAR` of any light attenuation model,
e.g. `TwoBandPhotosyntheticallyActiveRadiation(grid, PARFromShortwave(Qsw))`.

When coupled to NumericalEarth, `PARFromShortwave(grid)` builds the surface field for you and the
coupled model writes the ocean's penetrating shortwave into it each step.
"""
struct PARFromShortwave{SS, FT}
                       surface_shortwave :: SS
    photosynthetic_fraction_of_shortwave :: FT
end

Adapt.adapt_structure(to, ad::PARFromShortwave) =
    PARFromShortwave(adapt(to, ad.surface_shortwave),
                     adapt(to, ad.photosynthetic_fraction_of_shortwave))

PARFromShortwave(surface_shortwave;
                 photosynthetic_fraction_of_shortwave = 0.43) =
    PARFromShortwave(surface_shortwave, photosynthetic_fraction_of_shortwave)

@inline Oceananigans.BoundaryConditions.getbc(light::PARFromShortwave, i, j, args...) =
    @inbounds light.photosynthetic_fraction_of_shortwave * getbc(light.surface_shortwave, i, j, args...)

summary(::PARFromShortwave) = "PARFromShortwave"

function show(io::IO, light::PARFromShortwave)
    msg = "PARFromShortwave\n"
    msg *= "└── photosynthetic_fraction_of_shortwave: $(light.photosynthetic_fraction_of_shortwave)\n"

    print(io, msg)

    return nothing
end
