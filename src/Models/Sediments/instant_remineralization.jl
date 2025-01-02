using Oceananigans.Fields: tracernames

"""
    InstantRemineralisation

Hold the parameters and fields the simplest benthic boundary layer where
organic carbon is assumed to remineralise instantly with some portion 
becoming N, and a fraction being permanently buried.

Burial efficiency from [RemineralisationFraction](@citet).
"""
struct InstantRemineralisation{FT, ST, RR} <: AbstractContinuousFormSedimentBiogeochemistry
          burial_efficiency_constant1 :: FT
          burial_efficiency_constant2 :: FT
    burial_efficiency_half_saturation :: FT

                      sinking_tracers :: ST
            remineralisation_reciever :: RR

    function InstantRemineralisation(a::FT, b::FT, k::FT, 
                                     sinking_tracers::ST, 
                                     remineralisation_reciever::RR) where {FT, ST, RR}

        add_remineralisation_methods!(remineralisation_reciever)

        return new{FT, ST, RR}(a, b, k, sinking_tracers, remineralisation_reciever)
    end
end

"""
    InstantRemineralisationSediment(grid;
                                    burial_efficiency_constant1::FT = 0.013,
                                    burial_efficiency_constant2::FT = 0.53,
                                    burial_efficiency_half_saturation::FT = 7)

Return a single-layer instant remineralisaiton model for NPZD bgc models.

Example
=======

```@example
using OceanBioME, Oceananigans, OceanBioME.Sediments

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment_model = InstantRemineralisation(grid)
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
@inline sinking_fluxs(s::InstantRemineralisation) = s.sinking_tracers
@inline coupled_tracers(s::InstantRemineralisation) = tuple(s.remineralisation_reciever)

@inline function (s::InstantRemineralisation)(::Val{:storage}, x, y, t, args...)
    a = s.burial_efficiency_constant1
    b = s.burial_efficiency_constant2
    k = s.burial_efficiency_half_saturation

    flux = @inbounds sum(args[2:end])

    burial_efficiency = a + b * (flux / (k + flux)) ^ 2
    
    return burial_efficiency * flux
end

@inline function remineralisation(s::InstantRemineralisation, x, y, t, args...)
    a = s.burial_efficiency_constant1
    b = s.burial_efficiency_constant2
    k = s.burial_efficiency_half_saturation

    flux = @inbounds sum(args[2:end])

    burial_efficiency = a + b * (flux / (k + flux)) ^ 2

    return (1 - burial_efficiency) * flux
end

function add_remineralisation_methods!(remineralisation_reciever; fname = remineralisation)
    method = quote
        function (s::InstantRemineralisation)(::$(typeof(Val(remineralisation_reciever))), args...)
            return $(fname)(s, args...)
        end
    end

    eval(method)
end

summary(::InstantRemineralisation{FT}) where {FT} = string("Single-layer instant remineralisaiton ($FT)")
show(io::IO, model::InstantRemineralisation) = print(io, summary(model))
