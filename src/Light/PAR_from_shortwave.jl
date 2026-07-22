import Oceananigans.BoundaryConditions: getbc

struct PARFromShortwave{SN, FT}
                   surface_net_radiation :: SN
    photosynthetic_fraction_of_shortwave :: FT
end

Adapt.adapt_structure(to, ad::PARFromShortwave) =
    PARFromShortwave(adapt(to, ad.surface_net_radiation),
                                    nothing)

function PARFromShortwave(grid;
                          photosynthetic_fraction_of_shortwave = 0.43)

    surface_net_radiation = Field{Center, Center, Nothing}(grid)

    return PARFromShortwave(surface_net_radiation, photosynthetic_fraction_of_shortwave)
end

@inline getbc(light::PARFromShortwave, i, j, grid, clock, field) = 
    @inbounds light.surface_net_radiation[i, j, 1]
    