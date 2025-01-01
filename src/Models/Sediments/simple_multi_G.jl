import Base: show, summary

"""
    struct SimpleMultiG

Hold the parameters and fields for a simple "multi G" single-layer sediment model.
Based on the Level 3 model described by [Soetaert2000](@citet).
"""
struct SimpleMultiG{FT, P1, P2, P3, P4, F, TE, B} <: FlatSediment
             fast_decay_rate :: FT
             slow_decay_rate :: FT

               fast_redfield :: FT
               slow_redfield :: FT

               fast_fraction :: FT
               slow_fraction :: FT
          refactory_fraction :: FT

    nitrate_oxidation_params :: P1
      denitrification_params :: P2
               anoxic_params :: P3
            solid_dep_params :: P4

                      fields :: F
                  tendencies :: TE
              bottom_indices :: B

    SimpleMultiG(fast_decay_rate::FT, slow_decay_rate::FT,
                 fast_redfield::FT,   slow_redfield::FT,
                 fast_fraction::FT,   slow_fraction::FT, refactory_fraction::FT,
                 nitrate_oxidation_params::P1,
                 denitrification_params::P2,
                 anoxic_params::P3,
                 solid_dep_params::P4,
                 fields::F, tendencies::TE,
                 bottom_indices::B) where {FT, P1, P2, P3, P4, F, TE, B} =
        new{FT, P1, P2, P3, P4, F, TE, B}(fast_decay_rate, slow_decay_rate,
                                          fast_redfield,   slow_redfield,
                                          fast_fraction,   slow_fraction, refactory_fraction,
                                          nitrate_oxidation_params,
                                          denitrification_params,
                                          anoxic_params,
                                          solid_dep_params,
                                          fields, tendencies,
                                          bottom_indices)
end

"""
    SimpleMultiG(; grid
                   fast_decay_rate = 2/day,
                   slow_decay_rate = 0.2/day,
                   fast_redfield = 0.1509,
                   slow_redfield = 0.13,
                   fast_fraction = 0.74,
                   slow_fraction = 0.26,
                   refactory_fraction = 0.1,
                   nitrate_oxidation_params = on_architecture(architecture(grid), [- 1.9785, 0.2261, -0.0615, -0.0289, - 0.36109, - 0.0232]),
                   denitrification_params = on_architecture(architecture(grid), [- 3.0790, 1.7509, 0.0593, - 0.1923, 0.0604, 0.0662]),
                   anoxic_params = on_architecture(architecture(grid), [- 3.9476, 2.6269, - 0.2426, -1.3349, 0.1826, - 0.0143]),
                   solid_dep_params = on_architecture(architecture(grid), [0.233, 0.336, 982.0, - 1.548]))

Return a single-layer "multi G" sediment model (`SimpleMultiG`) on `grid`, where parameters
can be optionally specified.

The model is a single layer (i.e. does not include porous diffusion) model with three classes
of sediment organic matter which decay at three different rates (fast, slow, refactory).
The nitrification/denitrification/anoxic mineralisation fractions default to the parameterisation
of Soetaert et al. 2000; doi:[10.1016/S0012-8252(00)00004-0](https://doi.org/10.1016/S0012-8252(00)00004-0).

This model has not yet been validated or compared to observational data. The variety of degridation
processes is likely to be strongly dependent on oxygen availability (see
[https://bg.copernicus.org/articles/6/1273/2009/bg-6-1273-2009.pdf](https://bg.copernicus.org/articles/6/1273/2009/bg-6-1273-2009.pdf))
so it will therefore be important to also thoroughly validate the oxygen model (also currently limited).

Example
=======

```jldoctest simplemultig; filter = r".*@ OceanBioME.Models.Sediments.*"
julia> using OceanBioME, Oceananigans, OceanBioME.Sediments

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> sediment_model = SimpleMultiG(; grid)
┌ Warning: Sediment models are an experimental feature and have not yet been validated.
└ @ OceanBioME.Models.Sediments ~/Documents/Projects/OceanBioME.jl/src/Models/Sediments/simple_multi_G.jl:104
[ Info: This sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC.
Single-layer multi-G sediment model (Float64)
```
"""
function SimpleMultiG(; grid,
                        fast_decay_rate = 2/day,
                        slow_decay_rate = 0.2/day,
                        fast_redfield = 0.1509,
                        slow_redfield = 0.13,
                        fast_fraction = 0.74,
                        slow_fraction = 0.26,
                        refactory_fraction = 0.1,
                        nitrate_oxidation_params = on_architecture(architecture(grid), [- 1.9785, 0.2261, -0.0615, -0.0289, - 0.36109, - 0.0232]),
                        denitrification_params = on_architecture(architecture(grid), [- 3.0790, 1.7509, 0.0593, - 0.1923, 0.0604, 0.0662]),
                        anoxic_params = on_architecture(architecture(grid), [- 3.9476, 2.6269, - 0.2426, -1.3349, 0.1826, - 0.0143]),
                        solid_dep_params = on_architecture(architecture(grid), [0.233, 0.336, 982.0, - 1.548]))

    @warn "Sediment models are an experimental feature and have not yet been validated."
    @warn "Sediment models currently do not pass tests and are probably broken!"
    @info "The SimpleMultiG sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC."

    tracer_names = (:C_slow, :C_fast, :N_slow, :N_fast, :C_ref, :N_ref)

    # add field slicing back ( indices = (:, :, 1)) when output writer can cope
    fields = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names))
    tendencies = (Gⁿ = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names)),
                  G⁻ = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names)))

    bottom_indices = on_architecture(architecture(grid), calculate_bottom_indices(grid))

    return SimpleMultiG(fast_decay_rate, slow_decay_rate,
                        fast_redfield, slow_redfield,
                        fast_fraction, slow_fraction, refactory_fraction,
                        nitrate_oxidation_params,
                        denitrification_params,
                        anoxic_params,
                        solid_dep_params,
                        fields,
                        tendencies,
                        bottom_indices)
end

adapt_structure(to, sediment::SimpleMultiG) =
    SimpleMultiG(adapt(to, sediment.fast_decay_rate),
                 adapt(to, sediment.slow_decay_rate),
                 adapt(to, sediment.fast_redfield),
                 adapt(to, sediment.slow_redfield),
                 adapt(to, sediment.fast_fraction),
                 adapt(to, sediment.slow_fraction),
                 adapt(to, sediment.refactory_fraction),
                 adapt(to, sediment.nitrate_oxidation_params),
                 adapt(to, sediment.denitrification_params),
                 adapt(to, sediment.anoxic_params),
                 adapt(to, sediment.solid_dep_params),
                 adapt(to, sediment.fields),
                 nothing,
                 adapt(to, sediment.bottom_indices))

sediment_tracers(::SimpleMultiG) = (:C_slow, :C_fast, :C_ref, :N_slow, :N_fast, :N_ref)
sediment_fields(model::SimpleMultiG) = (C_slow = model.fields.C_slow,
                                        C_fast = model.fields.C_fast,
                                        N_slow = model.fields.N_slow,
                                        N_fast = model.fields.N_fast,
                                         C_ref = model.fields.C_ref,
                                         N_ref = model.fields.N_ref)

@inline required_tracers(::SimpleMultiG, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, sinking_tracers(bgc)...)]
@inline required_tendencies(::SimpleMultiG, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, :DIC)]

@inline bottom_index_array(sediment::SimpleMultiG) = sediment.bottom_indices

function _calculate_sediment_tendencies!(i, j, sediment::SimpleMultiG, bgc, grid, advection, tracers, tendencies, sediment_tendencies, t)
    k = bottom_index(i, j, sediment)
    depth = @inbounds -znodes(grid, Center(), Center(), Center())[k]

    Δz = zspacing(i, j, k, grid, Center(), Center(), Center())

    @inbounds begin
        carbon_deposition = carbon_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        nitrogen_deposition = nitrogen_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        # rates
        C_min_slow = sediment.fields.C_slow[i, j, 1] * sediment.slow_decay_rate
        C_min_fast = sediment.fields.C_fast[i, j, 1] * sediment.fast_decay_rate

        N_min_slow = sediment.fields.N_slow[i, j, 1] * sediment.slow_decay_rate
        N_min_fast = sediment.fields.N_fast[i, j, 1] * sediment.fast_decay_rate

        Cᵐⁱⁿ = C_min_slow + C_min_fast
        Nᵐⁱⁿ = N_min_slow + N_min_fast

        reactivity = Cᵐⁱⁿ * day / (sediment.fields.C_slow[i, j, 1] + sediment.fields.C_fast[i, j, 1])

        # sediment evolution
        sediment_tendencies.C_slow[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.slow_fraction * carbon_deposition - C_min_slow
        sediment_tendencies.C_fast[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.fast_fraction * carbon_deposition - C_min_fast
        sediment_tendencies.C_ref[i, j, 1] = sediment.refactory_fraction * carbon_deposition

        sediment_tendencies.N_slow[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.slow_fraction * nitrogen_deposition - N_min_slow
        sediment_tendencies.N_fast[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.fast_fraction * nitrogen_deposition - N_min_fast
        sediment_tendencies.N_ref[i, j, 1] = sediment.refactory_fraction * nitrogen_deposition

        # efflux/influx
        O₂  = tracers.O₂[i, j, k]
        NO₃ = tracers.NO₃[i, j, k]
        NH₄ = tracers.NH₄[i, j, k]

        A, B, C, D, E, F = sediment.nitrate_oxidation_params

        pₙᵢₜ = exp(A +
                  B * log(Cᵐⁱⁿ * day) * log(O₂) +
                  C * log(Cᵐⁱⁿ * day) ^ 2 +
                  D * log(reactivity) * log(NH₄) +
                  E * log(Cᵐⁱⁿ * day) +
                  F * log(Cᵐⁱⁿ * day) * log(NH₄)) / (Nᵐⁱⁿ * day)

        #=
        pᵈᵉⁿⁱᵗ = exp(sediment.denitrification_params.A +
                     sediment.denitrification_params.B * log(Cᵐⁱⁿ * day) +
                     sediment.denitrification_params.C * log(NO₃) ^ 2 +
                     sediment.denitrification_params.D * log(Cᵐⁱⁿ * day) ^ 2 +
                     sediment.denitrification_params.E * log(reactivity) ^ 2 +
                     sediment.denitrification_params.F * log(O₂) * log(reactivity)) / (Cᵐⁱⁿ * day)
        =#

        A, B, C, D, E, F = sediment.anoxic_params

        pₐₙₒₓ = exp(A +
                    B * log(Cᵐⁱⁿ * day) +
                    C * log(Cᵐⁱⁿ * day) ^ 2 +
                    D * log(reactivity) +
                    E * log(O₂) * log(reactivity) +
                    F * log(NO₃) ^ 2) / (Cᵐⁱⁿ * day)

        A, B, C, D = sediment.solid_dep_params
        pₛₒₗᵢ = A * (C * depth ^ D) ^ B

        tendencies.NH₄[i, j, k] += Nᵐⁱⁿ * (1 - pₙᵢₜ) / Δz
        tendencies.NO₃[i, j, k] += Nᵐⁱⁿ * pₙᵢₜ / Δz
        tendencies.DIC[i, j, k] += Cᵐⁱⁿ / Δz
        tendencies.O₂[i, j, k]  -= max(0, ((1 - pₐₙₒₓ * pₛₒₗᵢ) * Cᵐⁱⁿ + 2 * Nᵐⁱⁿ * pₙᵢₜ)/ Δz) # this seems dodge but this model doesn't cope with anoxia properly (I think)
    end
end

summary(::SimpleMultiG{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G sediment model ($FT)")
show(io::IO, model::SimpleMultiG) = print(io, summary(model))
