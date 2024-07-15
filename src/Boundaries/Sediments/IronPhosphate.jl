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

    tracer_names = (:O2, :NH4, :NO3, :NO2, :N2, :TPO4, :FeOHP, :FeII, :FeS2, :SO4, :TH2S, :CH4, :TCO2, :Gi)
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
                  
sediment_tracers(::IronPhosphate) = (:O2, :NH4, :NO3, :NO2, :N2, :TPO4, :FeOHP, :FeII, :FeS2, :SO4, :TH2S, :CH4, :TCO2, :Gi)
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

        carbon_deposition = carbon_flux(i, j, k, grid, advection, bgc, tracers) * Δz
                        
        nitrogen_deposition = nitrogen_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        phosphate_deposition = phosphate_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        iron_deposition = phosphate_flux(i, j, k, grid, advection, bgc, tracers) * Δz * 0.1 # molar ratio from Dale et al 2012, p. 633
        
        O2 = sediment.fields.O2[i, j, 1]
        NH4 = sediment.fields.NH4[i, j, 1]
        NO3 = sediment.fields.NO3[i, j, 1]
        NO2 = sediment.fields.NO2[i, j, 1]
        N2 = sediment.fields.N2[i, j, 1]
        TPO4 = sediment.fields.TPO4[i, j, 1]
        FeOHP = sediment.fields.FeOHP[i, j, 1]
        FeII = sediment.fields.FeII[i, j, 1]
        FeS2 = sediment.fields.FeS2[i, j, 1]
        SO4 = sediment.fields.SO4[i, j, 1]
        TH2S = sediment.fields.TH2S[i, j, 1]
        CH4 = sediment.fields.CH4[i, j, 1]
        TCO2 = sediment.fields.TCO2[i, j, 1]
        Gi = sediment.fields.Gi[i, j, 1]

        #####
        ##### RATES
        #####

        @inline KO2 = 1 # μM, Half–saturation constant for O2
        @inline KNO3 = 10 # μM, Half–saturation constant for NO3
        @inline KNO2 = 10 # μM, Half–saturation constant for NO2
        @inline KFe = 0.028 # wt-%, Half–saturation constant for Fe
        @inline KSO4 = 0.1 # μM, Half–saturation constant for SO4
        @inline KTPO4 = 10 # μM, Half–saturation constant for TPO4
        
        @inline kGi = 0.016 # day⁻¹, Rate constant for G0 degradation, Dale et al
        @inline fT = 1 # TODO temperature correction for rates
        @inline fox = 10 # Enhancement factor for POM degradation by O2

        fK_O2 = O2 / (O2 + KO2) # kinetic limiting term
        fK_NO3 = NO3 / (NO3 + KNO3) # kinetic limiting term
        fK_NO2 = NO2 / (NO2 + KNO2) # kinetic limiting term
        fK_Fe = FeOHP / (FeOHP + KFe) # kinetic limiting term
        fK_SO4 = SO4 / (SO4 + KSO4) # kinetic limiting term
        fK_TPO4 = TPO4 / (TPO4 + KTPO4) # kinetic limiting term

        RO2 = Gi * (fT * kGi * fox * fK_O2)
        RNO3 = Gi * (fT * kGi * fK_NO3 * (1 - fK_NO2) * (1 - fk_O2))
        RNO2 = Gi * (fT * kGi * fK_NO2 * (1 - fk_O2))
        RFe = Gi * (fT * kGi * fK_Fe * (1 - fK_NO3) * (1 - fK_NO2) * (1 - fk_NO2))
        RSO4 = Gi * (fT * kGi * fK_SO4 * (1 - fK_Fe) *(1 - fK_NO3) * (1 - fK_NO2) * (1 - fK_NO2))
        RCH4 = Gi * (fT * kGi * (1 - fK_SO4) * (1 - fK_Fe) *(1 - fK_NO3) * (1 - fK_NO3) * (1 - fK_NO2))

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

        RDNRA = TH2S * NO3 * fT * kDNRA
        Ramx = NH4 * NO2 * fT * kamx
        RNH4ox = NH4 * O2 * fT * kNH4ox
        RNO2ox = NO2 * O2 * fT * kNO2ox
        RAOM = CH4 * SO4 * fT * kAOM
        RH2Sox = TH2S * O2 * fT * kH2Sox
        RFe2ox = FeII * O2 * fT * kFe2ox
        RFeS2ox = FeS2 * O2 * fT * kFeS2ox
        RFeS2p = TH2S * FeII * fT * kFeS2p
        RFe3red = TH2S ^ 0.5 * FeOHP * fT * kFe3red * (2 / (O2 + 2))

        ratio_NC = 9.5/106
        ratio_PC = 1/106
        ratio_FeP = 0.1

        #####
        ##### sediment evolution
        #####

        sediment_tendencies.Gi = carbon_deposition - RO2 - RNO3 - RNO2 - RFe - RSO4 - RCH4

        tendencies.PO4[i, j, k] = 0.1 * fK_TPO4 * RFe2ox

    end
end

summary(::IronPhosphate{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G + Iron + Phosphate sediment model ($FT)")
show(io::IO, model::IronPhosphate) = print(io, summary(model))
