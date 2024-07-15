import Base: show, summary

"""
    struct IronPhosphate

Hold the parameters and fields for a single-layer sediment model including Iron and Phosphate.
Based on the Level 3 model described by [Soetaert2000](@citet) and the Iron and Phosphate
interactions described in A. Dale et al 2012 Biogeosciences.
"""
struct IronPhosphate{FT, P1, P2, P3, P4, F, TE, B} <: FlatSediment
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
    
    IronPhosphate(fast_decay_rate::FT, slow_decay_rate::FT,
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
    IronPhosphate(; grid
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

Return a single-layer "multi G" + Iron + Phosphate sediment model (`SimpleMultiG`) on `grid`, where parameters
can be optionally specified.

The model is a single layer (i.e. does not include porous diffusion) model with three classes
of sediment organic matter which decay at three different rates (fast, slow, refactory).
The nitrification/denitrification/anoxic mineralisation fractions default to the parameterisation
of Soetaert et al. 2000; doi:[10.1016/S0012-8252(00)00004-0](https://doi.org/10.1016/S0012-8252(00)00004-0).

Additionally, Iron and Phosphate reactions are computed to adjust O2, NO3 and POC fluxes

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
Single-layer multi-G + Iron + Phosphate sediment model (Float64)
```
"""
function IronPhosphate(; grid,
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
    @info "This sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC."

    tracer_names = (:O₂, :NH₄, :NO₃, :NO₂, :N₂, :TPO₄, :FeOHP, :Feᴵᴵ, :FeS₂, :SO₄, :TH₂S, :CH₄, :TCO₂, :Gi)
        #:C_slow, :C_fast, :N_slow, :N_fast, :C_ref, :N_ref, :Fe_III, :Fe_II, :PO4_dissolved, :P_org)

    # add field slicing back ( indices = (:, :, 1)) when output writer can cope
    fields = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names))
    tendencies = (Gⁿ = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names)),
                  G⁻ = NamedTuple{tracer_names}(Tuple(CenterField(grid) for tracer in tracer_names)))

    bottom_indices = on_architecture(architecture(grid), calculate_bottom_indices(grid))

    return IronPhosphate(fast_decay_rate, slow_decay_rate,
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

adapt_structure(to, sediment::IronPhosphate) = 
    IronPhosphate(adapt(to, sediment.fast_decay_rate),
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
                  
sediment_tracers(::IronPhosphate) = (:O₂, :NH₄, :NO₃, :NO₂, :N₂, :TPO₄, :FeOHP, :Feᴵᴵ, :FeS₂, :SO₄, :TH₂S, :CH₄, :TCO₂, :Gi)
sediment_fields(model::IronPhosphate) = (O2 = model.fields.O2,
                                         NH4 = model.fields.NH4,
                                         NO3 = model.fields.NO3,
                                         NO2 = model.fields.NO2,
                                         N2 = model.fields.N2,
                                         TPO4 = model.fields.TPO4,
                                         FeOHP = model.fields.FeOHP,
                                         FeII = model.fields.FeII,
                                         FeS2 = model.fields.FeS2,
                                         SO4 = model.fields.SO4,
                                         TH2S = model.fields.TH2S,
                                         CH4 = model.fields.CH4,
                                         TCO2 = model.fields.TCO2,
                                         Gi = model.fields.Gi)

@inline required_tracers(::IronPhosphate, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, :DIC, :POC, :Feᴾ, sinking_tracers(bgc)...)]
@inline required_tendencies(::IronPhosphate, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, :DIC, :POC, :Feᴾ)] # TODO add effluxes of PO4 and Fe once PISCES is working

@inline bottom_index_array(sediment::IronPhosphate) = sediment.bottom_indices

function _calculate_sediment_tendencies!(i, j, sediment::IronPhosphate, bgc, grid, advection, tracers, tendencies, sediment_tendencies, t)
    k = bottom_index(i, j, sediment)
    depth = @inbounds -znodes(grid, Center(), Center(), Center())[k]

    Δz = zspacing(i, j, k, grid, Center(), Center(), Center())

    @inbounds begin
        oxygen_deposition = oxygen_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        POC_deposition = poc_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        iron_deposition = iron_flux(i, j, k, grid, advection, bgc, tracers) * Δz
        
        #(:O₂, :NH₄, :NO₃, :NO₂, :N₂, :TPO₄, :FeOHP, :Feᴵᴵ, :FeS₂, :SO₄, :TH₂S, :CH₄, :TCO₂, :Gi)
        O₂ = sediment.fields.O₂[i, j, 1]
        NH₄ = sediment.fields.NH₄[i, j, 1]
        NO₃ = sediment.fields.NO₃[i, j, 1]
        NO₂ = sediment.fields.NO₂[i, j, 1]
        TPO₄ = sediment.fields.TPO₄[i, j, 1]
        FeOHP = sediment.fields.FeOHP[i, j, 1]
        Feᴵᴵ = sediment.fields.Feᴵᴵ[i, j, 1]
        FeS₂ = sediment.fields.FeS₂[i, j, 1]
        SO₄ = sediment.fields.SO₄[i, j, 1]
        TH₂S = sediment.fields.TH₂S[i, j, 1]
        CH₄ = sediment.fields.CH₄[i, j, 1]
        Gi = sediment.fields.Gi[i, j, 1]

        #####
        ##### RATES
        #####

        @inline K_O₂ = 1 # μM, Half–saturation constant for O2
        @inline K_NO₃ = 10 # μM, Half–saturation constant for NO3
        @inline K_NO₂ = 10 # μM, Half–saturation constant for NO2
        @inline K_Fe = 0.028 # wt-%, Half–saturation constant for Fe
        @inline K_SO₄ = 0.1 # μM, Half–saturation constant for SO4
        @inline K_TPO₄ = 10 # μM, Half–saturation constant for TPO4
        #(:O₂, :NH₄, :NO₃, :NO₂, :N₂, :TPO₄, :FeOHP, :Feᴵᴵ, :FeS₂, :SO₄, :TH₂S, :CH₄, :TCO₂, :Gi)
         
        @inline kGi = 0.016 # day⁻¹, Rate constant for G0 degradation, Dale et al
        @inline fT = 1 # TODO temperature correction for rates
        @inline fₒₓ = 10 # Enhancement factor for POM degradation by O2
        
        fₖ₋ₒ₂ = O₂ / (O₂ + K_O₂) # kinetic limiting term
        fₖ₋ₙₒ₃ = NO₃ / (NO₃ + K_NO₃) # kinetic limiting term
        fₖ₋ₙₒ₂ = NO₂ / (NO₂ + K_NO₂) # kinetic limiting term
        fK_Fe = FeOHP / (FeOHP + K_Fe) # kinetic limiting term
        fₖ₋ₛₒ₄ = SO₄ / (SO₄ + K_SO₄) # kinetic limiting term
        fₖ₋ₜₚₒ₄ = TPO₄ / (TPO₄ + K_TPO₄) # kinetic limiting term

        RO₂ = Gi * (fT * kGi * fₒₓ * fₖ₋ₒ₂)
        RNO₃ = Gi * (fT * kGi * fₖ₋ₙₒ₃ * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₒ₂))
        RNO₂ = Gi * (fT * kGi * fₖ₋ₙₒ₂ * (1 - fₖ₋ₒ₂))
        RFe = Gi * (fT * kGi * fK_Fe * (1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₙₒ₂))
        RSO₄ = Gi * (fT * kGi * fₖ₋ₛₒ₄ * (1 - fK_Fe) *(1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₙₒ₂))
        RCH₄ = Gi * (fT * kGi * (1 - fₖ₋ₛₒ₄) * (1 - fK_Fe) *(1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂))

        @inline kDNRA = 2.7e5 # M-1 day-1
        @inline kamx = 2.7e4 # M-1 day-1
        @inline kNH4ox = 2.7e4# M-1 day-1
        @inline kNO2ox = 2.7e4 # M-1 day-1
        @inline kAOM = 0.27 # day-1
        @inline kH2Sox = 2.7e4 # M-1 day-1
        @inline kFe2ox = 2.7e5 # M-1 day-1
        @inline kFeS2ox = 2.7e3 # M-1 day-1
        @inline kFeS2p = 2.7e4 # M-1 day-1
        @inline kFe3red = 0.82 # cm1.5 mmol−0.5 day−1

        RDNRA = TH₂S * NO₃ * fT * kDNRA
        Ramx = NH₄ * NO₂ * fT * kamx
        RNH4ox = NH₄ * O₂ * fT * kNH4ox
        RNO2ox = NO₂ * O₂ * fT * kNO2ox
        RAOM = CH₄ * SO₄ * fT * kAOM
        RH2Sox = TH₂S * O₂ * fT * kH2Sox
        RFe2ox = Feᴵᴵ * O₂ * fT * kFe2ox
        RFeS2ox = FeS₂ * O₂ * fT * kFeS2ox
        RFeS2p = TH₂S * Feᴵᴵ * fT * kFeS2p # H2 should be produced here but I assume it dissociates...
        RFe3red = TH₂S ^ 0.5 * FeOHP * fT * kFe3red * (2 / (O₂ + 2))

        @inline ratio_NC = 9.5/106
        @inline ratio_PC = 1/106
        @inline ratio_FeP = 0.1

        #####
        ##### sediment evolution
        #####

        sediment_tendencies.Gi = POC_deposition - RO₂ - RNO₃ - RNO₂ - RFe - RSO₄ - RCH₄
        sediment_tendencies.O₂ = oxygen_deposition - RO₂ - 1.5 * RNH4ox - 0.5 * RNO2ox - 2 * RH2Sox - 0.25 * RFe2ox - 3.5 * RFeS2ox
        sediment_tendencies.TCO₂ = RO₂ + RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄ + RAOM
        sediment_tendencies.NH₄ = ratio_NC * (RO₂ + RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄) + RDNRA - Ramx - RNH4ox
        sediment_tendencies.NO₃ = -2 * RNO₃ - RDNRA + RNO2ox
        sediment_tendencies.NO₂ = 2 * RNO₃ - 1.33 * RNO₂ - Ramx + RNH4ox - RNO2ox
        sediment_tendencies.N₂ = Ramx + 0.66 * RNO₂
        sediment_tendencies.TPO₄ = ratio_PC * (RO₂ + RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄) + RFe * ratio_FeP - ratio_FeP * RFe2ox * fₖ₋ₜₚₒ₄ + ratio_FeP * RFe3red
        sediment_tendencies.FeOHP = iron_deposition - 4 * RFe + RFe2ox - RFe3red
        sediment_tendencies.Feᴵᴵ = 4 * RFe - RFe2ox + RFeS2ox - RFeS2p + RFe3red
        sediment_tendencies.FeS₂ = RFeS2p - RFeS2ox
        sediment_tendencies.SO₄ = -0.5 * RSO₄ + RDNRA - RAOM + RH2Sox + 2 * RFeS2ox + (RFe3red / 8)
        sediment_tendencies.TH₂S = 0.5 * RSO₄ - RDNRA + RAOM - RH2Sox  - 2 * RFeS2p - (RFe3red / 8)
        sediment_tendencies.CH₄ = 0.5 * RCH₄ - RAOM
    end
end

summary(::IronPhosphate{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G + Iron + Phosphate sediment model ($FT)")
show(io::IO, model::IronPhosphate) = print(io, summary(model))
