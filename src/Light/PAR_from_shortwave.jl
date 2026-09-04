import Oceananigans.BoundaryConditions: getbc

struct PARFromShortwave{SS, FT}
                       surface_shortwave :: SS
    photosynthetic_fraction_of_shortwave :: FT
end

Adapt.adapt_structure(to, ad::PARFromShortwave) =
    PARFromShortwave(adapt(to, ad.surface_shortwave),
                     adapt(to, ad.surface_shortwave))

PARFromShortwave(surface_shortwave;
                 photosynthetic_fraction_of_shortwave = 0.43) = 
    PARFromShortwave(surface_shortwave, photosynthetic_fraction_of_shortwave)

@inline Oceananigans.BoundaryConditions.getbc(light::PARFromShortwave, i, j, args...) = 
    @inbounds light.photosynthetic_fraction_of_shortwave * getbc(light.surface_shortwave, i, j, args...)

Base.summary(io, light::PARFromShortwave) = "PARFromShortwave"

function Base.show(io::IO, light::PARFromShortwave) 
    msg = "PARFromShortwave\n"
    msg *= "└── photosynthetic_fraction_of_shortwave: $(light.photosynthetic_fraction_of_shortwave)\n"

    print(io, msg)

    return nothing
end