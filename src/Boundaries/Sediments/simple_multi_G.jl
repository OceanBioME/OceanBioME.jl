import Base: show, summary

"""
    struct SimpleMultiG

Hold the parameters and fields for a simple "multi G" single-layer sediment model.
Based on the Level 3 model described by [Soetaert2000](@citet).
"""
struct SimpleMultiG{FT, P1, P2, P3, P4, F, TE} <: FlatSediment
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
end

"""
    SimpleMultiG(; grid
                   fast_decay_rate::FT = 2/day,
                   slow_decay_rate::FT = 0.2/day,
                   fast_redfield::FT = 0.1509,
                   slow_redfield::FT = 0.13,
                   fast_fraction::FT = 0.74,
                   slow_fraction::FT = 0.26,
                   refactory_fraction::FT = 0.1,
                   nitrate_oxidation_params::P1 = (A = - 1.9785, B = 0.2261, C = -0.0615, D = -0.0289, E = - 0.36109, F = - 0.0232),
                   denitrification_params::P2 = (A = - 3.0790, B = 1.7509, C = 0.0593, D = - 0.1923, E = 0.0604, F = 0.0662),
                   anoxic_params::P3 = (A = - 3.9476, B = 2.6269, C = - 0.2426, D = -1.3349, E = 0.1826, F = - 0.0143),
                   depth = abs(znodes(grid, Face())[1]),
                   solid_dep_params::P4 = (A = 0.233, B = 0.336, C = 982, D = - 1.548, depth = depth))

Return a single-layer "multi G" sediment model (`SimpleMultiG`) on `grid`, where parameters
can be optionally specified.

The model is a single layer (i.e. does not include porous diffusion) model with three classes
of sediment organic matter which decay at three different rates (fast, slow, refactory).
The nitrifcation/denitrifcation/anoxic mineralisation fractions default to the parameterisation
of Soetaert et al. 2000; doi:[10.1016/S0012-8252(00)00004-0](https://doi.org/10.1016/S0012-8252(00)00004-0).

This model has not yet been validated or compared to observational data. The variety of degridation
processes is likely to be strongly dependent on oxygen availability (see
[https://bg.copernicus.org/articles/6/1273/2009/bg-6-1273-2009.pdf](https://bg.copernicus.org/articles/6/1273/2009/bg-6-1273-2009.pdf))
so it will therefore be important to also thoroughly validate the oxygen model (also currently limited).

Example
=======

```jldoctest simplemultig; filter = r".*@ OceanBioME.Boundaries.Sediments.*"
julia> using OceanBioME, Oceananigans, OceanBioME.Sediments

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> sediment_model = SimpleMultiG(; grid)
┌ Warning: Sediment models are an experimental feature and have not yet been validated.
└ @ OceanBioME.Boundaries.Sediments ~/OceanBioME.jl/src/Boundaries/Sediments/simple_multi_G.jl:102
[ Info: This sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC.
Single-layer multi-G sediment model (Float64)
```
"""
function SimpleMultiG(; grid,
                       fast_decay_rate::FT = 2/day,
                       slow_decay_rate::FT = 0.2/day,
                       fast_redfield::FT = 0.1509,
                       slow_redfield::FT = 0.13,
                       fast_fraction::FT = 0.74,
                       slow_fraction::FT = 0.26,
                       refactory_fraction::FT = 0.1,
                       nitrate_oxidation_params::P1 = (A = - 1.9785,
                                                       B = 0.2261,
                                                       C = -0.0615,
                                                       D = -0.0289,
                                                       E = - 0.36109,
                                                       F = - 0.0232),
                       denitrification_params::P2 = (A = - 3.0790,
                                                     B = 1.7509,
                                                     C = 0.0593,
                                                     D = - 0.1923,
                                                     E = 0.0604,
                                                     F = 0.0662),
                       anoxic_params::P3 = (A  = - 3.9476,
                                            B = 2.6269,
                                            C = - 0.2426,
                                            D = -1.3349,
                                            E = 0.1826,
                                            F = - 0.0143),
                       depth = abs(znodes(grid, Face())[1]),
                       solid_dep_params::P4 = (A = 0.233, 
                                               B = 0.336, 
                                               C = 982, 
                                               D = - 1.548,
                                               depth = depth)) where {FT, P1, P2, P3, P4}

    @warn "Sediment models are an experimental feature and have not yet been validated."
    @info "This sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC."

    tracer_names = (:C_slow, :C_fast, :N_slow, :N_fast, :C_ref, :N_ref)

    # add field slicing back ( indices = (:, :, 1)) when output writer can cope
    fields = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names))
    tendencies = (Gⁿ = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names)),
                  G⁻ = NamedTuple{tracer_names}(Tuple(CenterField(grid;) for tracer in tracer_names)))

    F = typeof(fields)
    TE = typeof(tendencies)

    return SimpleMultiG{FT, P1, P2, P3, P4, F, TE}(fast_decay_rate, slow_decay_rate,
                                                   fast_redfield, slow_redfield,
                                                   fast_fraction, slow_fraction, refactory_fraction,
                                                   nitrate_oxidation_params,
                                                   denitrification_params,
                                                   anoxic_params,
                                                   solid_dep_params,
                                                   fields,
                                                   tendencies)
end

adapt_structure(to, sediment::SimpleMultiG) = 
    SimpleMultiG(sediment.fast_decay_rate,
                 sediment.slow_decay_rate,
                 sediment.fast_redfield,
                 sediment.slow_redfield,
                 sediment.fast_fraction,
                 sediment.slow_fraction,
                 sediment.refactory_fraction,
                 sediment.nitrate_oxidation_params,
                 sediment.denitrification_params,
                 sediment.anoxic_params,
                 sediment.solid_dep_params,
                 adapt(to, sediment.fields),
                 adapt(to, sediment.tendencies))
                  
sediment_tracers(::SimpleMultiG) = (:C_slow, :C_fast, :C_ref, :N_slow, :N_fast, :N_ref)
sediment_fields(model::SimpleMultiG) = (C_slow = model.fields.C_slow,
                                        C_fast = model.fields.C_fast,
                                        N_slow = model.fields.N_slow,
                                        N_fast = model.fields.N_fast,
                                         C_ref = model.fields.C_ref,
                                         N_ref = model.fields.N_ref)

@kernel function _calculate_tendencies!(sediment::SimpleMultiG, bgc, grid, advection, tracers, timestepper)
    i, j = @index(Global, NTuple)

    Δz = zspacing(i, j, 1, grid, Center(), Center(), Center())

    @inbounds begin

        carbon_deposition = carbon_flux(grid, advection, bgc, tracers, i, j) * Δz
                        
        nitrogen_deposition = nitrogen_flux(grid, advection, bgc, tracers, i, j) * Δz

        # rates
        C_min_slow = sediment.fields.C_slow[i, j, 1] * sediment.slow_decay_rate
        C_min_fast = sediment.fields.C_fast[i, j, 1] * sediment.fast_decay_rate

        N_min_slow = sediment.fields.N_slow[i, j, 1] * sediment.slow_decay_rate
        N_min_fast = sediment.fields.N_fast[i, j, 1] * sediment.fast_decay_rate

        Cᵐⁱⁿ = C_min_slow + C_min_fast
        Nᵐⁱⁿ = N_min_slow + N_min_fast
        
        k = Cᵐⁱⁿ * day / (sediment.fields.C_slow[i, j, 1] + sediment.fields.C_fast[i, j, 1])

        # sediment evolution
        sediment.tendencies.Gⁿ.C_slow[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.slow_fraction * carbon_deposition - C_min_slow
        sediment.tendencies.Gⁿ.C_fast[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.fast_fraction * carbon_deposition - C_min_fast
        sediment.tendencies.Gⁿ.C_ref[i, j, 1] = sediment.refactory_fraction * carbon_deposition

        sediment.tendencies.Gⁿ.N_slow[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.slow_fraction * nitrogen_deposition - N_min_slow
        sediment.tendencies.Gⁿ.N_fast[i, j, 1] = (1 - sediment.refactory_fraction) * sediment.fast_fraction * nitrogen_deposition - N_min_fast
        sediment.tendencies.Gⁿ.N_ref[i, j, 1] = sediment.refactory_fraction * nitrogen_deposition

        # efflux/influx
        O₂  = tracers.O₂[i, j, 1]
        NO₃ = tracers.NO₃[i, j, 1]
        NH₄ = tracers.NH₄[i, j, 1]

        pₙᵢₜ = exp(sediment.nitrate_oxidation_params.A +
                   sediment.nitrate_oxidation_params.B * log(Cᵐⁱⁿ * day) * log(O₂) +
                   sediment.nitrate_oxidation_params.C * log(Cᵐⁱⁿ * day) ^ 2 +
                   sediment.nitrate_oxidation_params.D * log(k) * log(NH₄) +
                   sediment.nitrate_oxidation_params.E * log(Cᵐⁱⁿ * day) +
                   sediment.nitrate_oxidation_params.F * log(Cᵐⁱⁿ * day) * log(NH₄)) / (Nᵐⁱⁿ * day)

        pᵈᵉⁿⁱᵗ = exp(sediment.denitrification_params.A +
                     sediment.denitrification_params.B * log(Cᵐⁱⁿ * day) +
                     sediment.denitrification_params.C * log(NO₃) ^ 2 +
                     sediment.denitrification_params.D * log(Cᵐⁱⁿ * day) ^ 2 +
                     sediment.denitrification_params.E * log(k) ^ 2 +
                     sediment.denitrification_params.F * log(O₂) * log(k)) / (Cᵐⁱⁿ * day)

        pₐₙₒₓ = exp(sediment.anoxic_params.A +
                    sediment.anoxic_params.B * log(Cᵐⁱⁿ * day) +
                    sediment.anoxic_params.C * log(Cᵐⁱⁿ * day) ^ 2 +
                    sediment.anoxic_params.D * log(k) +
                    sediment.anoxic_params.E * log(O₂) * log(k) +
                    sediment.anoxic_params.F * log(NO₃) ^ 2) / (Cᵐⁱⁿ * day)

        if isnan(pₐₙₒₓ)
            println("$(Cᵐⁱⁿ), $(k), $(O₂), $(NO₃)")
            error("Sediment anoxia has caused model failure")
        end

        pₛₒₗᵢ = sediment.solid_dep_params.A * (sediment.solid_dep_params.C * sediment.solid_dep_params.depth ^ sediment.solid_dep_params.D) ^ sediment.solid_dep_params.B

        Δz = grid.Δzᵃᵃᶜ[1]

        timestepper.Gⁿ.NH₄[i, j, 1] += Nᵐⁱⁿ * (1 - pₙᵢₜ) / Δz
        timestepper.Gⁿ.NO₃[i, j, 1] += Nᵐⁱⁿ * pₙᵢₜ / Δz
        timestepper.Gⁿ.DIC[i, j, 1] += Cᵐⁱⁿ / Δz
        timestepper.Gⁿ.O₂[i, j, 1]  -= max(0, (1 - pᵈᵉⁿⁱᵗ - pₐₙₒₓ * pₛₒₗᵢ) * Cᵐⁱⁿ / Δz) # this seems dodge but this model doesn't cope with anoxia properly
    end
end

summary(::SimpleMultiG{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G sediment model ($FT)")
show(io::IO, model::SimpleMultiG) = print(io, summary(model))
