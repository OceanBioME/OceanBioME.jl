using Oceananigans.Fields: tracernames

"""
    InstantRemineralisation

Hold the parameters and fields a simple sediment model where sinking organic
carbon is "instantly remineralised" and either returned to the domain as 
`remineralisation_reciever` (typically NH₄), or permanently stored in the sediment.

The "burial efficiency" (the fraction permanently stored) is from 
[RemineralisationFraction](@citet), and varies with the sinking flux.
"""
struct InstantRemineralisation{FT, ST, RR} <: AbstractContinuousFormSedimentBiogeochemistry
          burial_efficiency_constant1 :: FT
          burial_efficiency_constant2 :: FT
    burial_efficiency_half_saturation :: FT

                      sinking_tracers :: ST
end

# `RR` is the name of the tracer that receives the remineralised flux. It is a type parameter
# rather than a field so that the tendency method for that tracer (below) can dispatch on it.
InstantRemineralisation(a::FT, b::FT, k::FT, sinking_tracers::ST, remineralisation_reciever::Symbol) where {FT, ST} =
    InstantRemineralisation{FT, ST, remineralisation_reciever}(a, b, k, sinking_tracers)

@inline remineralisation_reciever(::InstantRemineralisation{<:Any, <:Any, RR}) where RR = RR

function Adapt.adapt_structure(to, ir::InstantRemineralisation{<:Any, <:Any, RR}) where RR
    a = adapt(to, ir.burial_efficiency_constant1)
    b = adapt(to, ir.burial_efficiency_constant2)
    k = adapt(to, ir.burial_efficiency_half_saturation)

    return InstantRemineralisation{typeof(a), Nothing, RR}(a, b, k, nothing)
end

"""
    InstantRemineralisationSediment(grid;
                                    sinking_tracers = (:P, :D), 
                                    remineralisation_reciever = :N,
                                    burial_efficiency_constant1 = 0.013,
                                    burial_efficiency_constant2 = 0.53,
                                    burial_efficiency_half_saturation = 7.0 / 6.56,
                                    kwargs...)

Return a single-layer instant remineralisation sediment model where the `sinking_tracers`
are instantly remineralised and returned to `remineralisation_reciever` with a small
fraction permanently buried with efficiency:

e = a + b * f / (k + f)²

where `a` is `burial_efficiency_constant1`, `b` is `burial_efficiency_constant2`, and 
`k` is the `burial_efficiency_half_saturation`.

`kwargs...` are `BiogeochemicalSediment` key word arguments.

Example
=======

```@example
using OceanBioME, Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment_model = InstantRemineralisationSediment(grid)

biogeochemistry = NPZD(; grid, sediment_model)
```

```@example
using OceanBioME, Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment_model = InstantRemineralisationSediment(grid; 
                                                 sinking_tracers = (:sPOM, :bPOM),
                                                 remineralisation_reciever = :NH₄)

biogeochemistry = LOBSTER(grid; sediment_model)
```
"""
InstantRemineralisationSediment(grid;
                                sinking_tracers = (:P, :D), 
                                remineralisation_reciever = :N,
                                burial_efficiency_constant1 = 0.013,
                                burial_efficiency_constant2 = 0.53,
                                burial_efficiency_half_saturation = 7.0 / 6.56,
                                kwargs...) =
    BiogeochemicalSediment(grid, 
                           InstantRemineralisation(burial_efficiency_constant1, 
                                                   burial_efficiency_constant2, 
                                                   burial_efficiency_half_saturation,
                                                   tracernames(sinking_tracers),
                                                   remineralisation_reciever);
                           kwargs...)
                  
@inline required_sediment_fields(::InstantRemineralisation) = (:storage, )
@inline required_tracers(::InstantRemineralisation) = tuple()
@inline sinking_fluxes(s::InstantRemineralisation) = s.sinking_tracers
@inline coupled_tracers(s::InstantRemineralisation) = tuple(remineralisation_reciever(s))

@inline function (s::InstantRemineralisation)(::Val{:storage}, x, y, t, storage, fluxs...)
    a = s.burial_efficiency_constant1
    b = s.burial_efficiency_constant2
    k = s.burial_efficiency_half_saturation

    flux = sum(fluxs)

    burial_efficiency = a + b * (flux / (k + flux)) ^ 2
    
    return burial_efficiency * flux
end

@inline function remineralisation(s::InstantRemineralisation, x, y, t, storage, fluxs...)
    a = s.burial_efficiency_constant1
    b = s.burial_efficiency_constant2
    k = s.burial_efficiency_half_saturation

    flux = sum(fluxs)

    burial_efficiency = a + b * (flux / (k + flux)) ^ 2

    return (1 - burial_efficiency) * flux
end

# the receiver's tendency is the remineralised flux
@inline (s::InstantRemineralisation{<:Any, <:Any, RR})(::Val{RR}, x, y, t, storage, fluxs...) where RR =
    remineralisation(s, x, y, t, storage, fluxs...)

summary(::InstantRemineralisation{FT}) where {FT} = string("Single-layer instant remineralisation ($FT)")
show(io::IO, model::InstantRemineralisation) = print(io, summary(model))
