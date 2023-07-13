"""
    struct InstantRemineralisation

Hold the parameters and fields the simplest benthic boundary layer where
organic carbon is assumed to remineralise instantly with some portion 
becoming N, and a fraction being perminantly burried.

Burial efficiency by 10.1029/2006GB002907, 2007
"""
struct InstantRemineralisation{FT, F, TE} <: FlatSediment
          burial_efficiency_constant1 :: FT
          burial_efficiency_constant2 :: FT
    burial_efficiency_half_saturaiton :: FT

                               fields :: F
                           tendencies :: TE

    function InstantRemineralisation(burial_efficiency_constant1::FT,
                                     burial_efficiency_constant2::FT,
                                     burial_efficiency_half_saturaiton::FT,
                                     fields::F,
                                     tendencies::TE) where {FT, F, TE}

        return new{FT, F, TE}(burial_efficiency_constant1,
                              burial_efficiency_constant2,
                              burial_efficiency_half_saturaiton,
                              fields,
                              tendencies)
    end
end

"""
    InstantRemineralisation(; grid,
        burial_efficiency_constant1::FT = 0.013,
        burial_efficiency_constant2::FT = 0.53,
        burial_efficiency_half_saturaiton::FT = 7) where FT

Returns a single-layer instant remineralisaiton model for NPZD bgc models.

Example
=======

```@example
using OceanBioME, Oceananigans, OceanBioME.Sediments

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment_model = InstantRemineralisation(; grid)
```
"""
function InstantRemineralisation(; grid,
            burial_efficiency_constant1 = 0.013,
            burial_efficiency_constant2 = 0.53,
            burial_efficiency_half_saturaiton = 7.0)

    @warn "Sediment models are an experimental feature and have not yet been validated"

    tracer_names = (:N_storage, )

    # add field slicing back ( indices = (:, :, 1)) when output writer can cope
    fields = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names))
    tendencies = (Gⁿ = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names)),
                  G⁻ = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names)))

    return InstantRemineralisation(burial_efficiency_constant1,
                                   burial_efficiency_constant2,
                                   burial_efficiency_half_saturaiton,
                                   fields,
                                   tendencies)
end

adapt_structure(to, sediment::InstantRemineralisation) = 
    SimpleMultiG(sediment.burial_efficiency_constant1,
                 sediment.burial_efficiency_constant2,
                 sediment.burial_efficiency_half_saturaiton,
                 adapt(to, sediment.fields),
                 adapt(to, sediment.tendencies))
                  
sediment_tracers(::InstantRemineralisation) = (:N_storage, )
sediment_fields(model::InstantRemineralisation) = (N_storage = model.fields.N_storage, )

@kernel function _calculate_tendencies!(sediment::InstantRemineralisation, bgc, grid, advection, tracers, timestepper)
    i, j = @index(Global, NTuple)

    @inbounds begin                        
        Δz = zspacing(i, j, 1, grid, Center(), Center(), Center())

        flux = nitrogen_flux(grid, advection, bgc, tracers, i, j) * Δz

        burial_efficiency = sediment.burial_efficiency_constant1 + sediment.burial_efficiency_constant2 * ((flux / 6.56) / (7 + flux / 6.56)) ^ 2

        # sediment evolution
        sediment.tendencies.Gⁿ.N_storage[i, j, 1] = burial_efficiency * flux

        remineralizaiton_reciever(bgc, timestepper.Gⁿ)[i, j, 1] += flux * (1 - burial_efficiency) / Δz
    end
end