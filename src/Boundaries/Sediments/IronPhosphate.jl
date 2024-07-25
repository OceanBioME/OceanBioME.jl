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

julia> sediment_model = IronPhosphate(; grid)
┌ Warning: Sediment models are an experimental feature and have not yet been validated.
└ @ OceanBioME.Boundaries.Sediments ~/OceanBioME.jl/src/Boundaries/Sediments/simple_multi_G.jl:102
[ Info: This sediment model is currently only compatible with models providing lots of different things.
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
sediment_fields(model::IronPhosphate) = (O₂ = model.fields.O₂,
                                        NH₄ = model.fields.NH₄,
                                        NO₃ = model.fields.NO₃,
                                        NO₂ = model.fields.NO₂,
                                        TPO₄ = model.fields.TPO₄,
                                        FeOHP = model.fields.FeOHP,
                                        Feᴵᴵ = model.fields.Feᴵᴵ,
                                        FeS₂ = model.fields.FeS₂,
                                        SO₄ = model.fields.SO₄,
                                        TH₂S = model.fields.TH₂S,
                                        CH₄ = model.fields.CH₄,
                                        Gi = model.fields.Gi)

@inline required_tracers(::IronPhosphate, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, :DIC)] #, sinking_tracers(bgc)...)]
@inline required_tendencies(::IronPhosphate, bgc, tracers) = tracers[(:NO₃, :NH₄, :O₂, :DIC)] # TODO add effluxes of PO4 and Fe once PISCES is working

@inline bottom_index_array(sediment::IronPhosphate) = sediment.bottom_indices

function _calculate_sediment_tendencies!(i, j, sediment::IronPhosphate, bgc, grid, advection, tracers, tendencies, sediment_tendencies, t)
    k = bottom_index(i, j, sediment)
    depth = @inbounds -znodes(grid, Center(), Center(), Center())[k]

    Δz = zspacing(i, j, k, grid, Center(), Center(), Center())

    @inbounds begin
        oxygen_deposition = 0#oxygen_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        POC_deposition = 5e-9 #carbon_flux(i, j, k, grid, advection, bgc, tracers) * Δz

        iron_deposition = 0#5e-10 #carbon_flux(i, j, k, grid, advection, bgc, tracers) * Δz * (1/106) * 0.1
        
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

        #println(FeOHP + FeS₂ + Feᴵᴵ) # iron is conserved
        #println(2*FeS₂+SO₄+TH₂S) # sulfur is conserved
        #saprintln(TPO₄)
        #####
        ##### RATES
        #####

        @inline K_O₂ = 1e-6 # μM, Half–saturation constant for O2
        @inline K_NO₃ = 10e-6 # μM, Half–saturation constant for NO3
        @inline K_NO₂ = 10e-6 # μM, Half–saturation constant for NO2
        @inline K_Fe = 0.028 # wt-%, Half–saturation constant for Fe
        @inline K_SO₄ = 0.1e-6 # μM, Half–saturation constant for SO4
        @inline K_TPO₄ = 10e-6 # μM, Half–saturation constant for TPO4
        #(:O₂, :NH₄, :NO₃, :NO₂, :N₂, :TPO₄, :FeOHP, :Feᴵᴵ, :FeS₂, :SO₄, :TH₂S, :CH₄, :TCO₂, :Gi)
         
        @inline kGi = 0.016 # day⁻¹, Rate constant for G0 degradation, Dale et al
        @inline fT = 1 # TODO temperature correction for rates
        @inline fₒₓ = 10 # Enhancement factor for POM degradation by O2
        @inline per_day_to_per_seconds = 1 / (60 * 60 * 24) 
        
        fₖ₋ₒ₂ = O₂ / (O₂ + K_O₂) # kinetic limiting term
        fₖ₋ₙₒ₃ = NO₃ / (NO₃ + K_NO₃) # kinetic limiting term
        fₖ₋ₙₒ₂ = NO₂ / (NO₂ + K_NO₂) # kinetic limiting term
        fK_Fe = FeOHP / (FeOHP + K_Fe) # kinetic limiting term
        fₖ₋ₛₒ₄ = SO₄ / (SO₄ + K_SO₄) # kinetic limiting term
        fₖ₋ₜₚₒ₄ = TPO₄ / (TPO₄ + K_TPO₄) # kinetic limiting term

        RO₂ = max(0, Gi * (fT * kGi * fₒₓ * fₖ₋ₒ₂) * per_day_to_per_seconds)
        RNO₃ = max(0, Gi * (fT * kGi * fₖ₋ₙₒ₃ * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₒ₂)) * per_day_to_per_seconds)
        RNO₂ = max(0, Gi * (fT * kGi * fₖ₋ₙₒ₂ * (1 - fₖ₋ₒ₂)) * per_day_to_per_seconds)
        RFe = max(0, Gi * (fT * kGi * fK_Fe * (1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₙₒ₂)) * per_day_to_per_seconds)
        RSO₄ = max(0, Gi * (fT * kGi * fₖ₋ₛₒ₄ * (1 - fK_Fe) *(1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₙₒ₂)) * per_day_to_per_seconds)
        RCH₄ = max(0, Gi * (fT * kGi * (1 - fₖ₋ₛₒ₄) * (1 - fK_Fe) *(1 - fₖ₋ₙₒ₃) * (1 - fₖ₋ₙₒ₂) * (1 - fₖ₋ₙₒ₂)) * per_day_to_per_seconds)

        

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

        R_DNRA = max(0, TH₂S * NO₃ * fT * kDNRA * per_day_to_per_seconds)
        R_amx = max(0, NH₄ * NO₂ * fT * kamx * per_day_to_per_seconds)
        R_NH4ox = max(0, NH₄ * O₂ * fT * kNH4ox * per_day_to_per_seconds)
        R_NO2ox = max(0, NO₂ * O₂ * fT * kNO2ox * per_day_to_per_seconds)
        R_AOM = max(0, CH₄ * SO₄ * fT * kAOM * per_day_to_per_seconds)
        R_H2Sox = max(0, TH₂S * O₂ * fT * kH2Sox * per_day_to_per_seconds)
        R_Fe2ox = max(0, Feᴵᴵ * O₂ * fT * kFe2ox * per_day_to_per_seconds)
        R_FeS2ox = max(0, FeS₂ * O₂ * fT * kFeS2ox * per_day_to_per_seconds)
        R_FeS2p = max(0, TH₂S * Feᴵᴵ * fT * kFeS2p * per_day_to_per_seconds) # H2 should be produced here but I assume it dissociates...
        R_Fe3red = max(0, max(0, TH₂S) ^ 0.5 * FeOHP * fT * kFe3red * (2 / (O₂ + 2)) * per_day_to_per_seconds)

        #println(R_DNRA, "  ", R_amx, "  ", R_NH4ox, "  ", R_NO2ox, "  ", R_AOM, "  ", R_H2Sox, "  ", R_Fe2ox, "  ", R_FeS2ox, "  ", R_FeS2p, "  ", R_Fe3red, "  ")
        if isnan(TPO₄)
            sleep(10)
        end
        println(O₂)
        @inline ratio_NC = 9.5/106
        @inline ratio_PC = 1/106
        @inline ratio_FeP = 0.1

        #####
        ##### sediment evolution
        #####

        sediment_tendencies.Gi[i, j, 1] = POC_deposition - RO₂ - RNO₃ - RNO₂ - RFe - RSO₄ - RCH₄
        sediment_tendencies.O₂[i, j, 1] = oxygen_deposition - RO₂ - 1.5 * R_NH4ox - 0.5 * R_NO2ox - 2 * R_H2Sox - 0.25 * R_Fe2ox - 3.5 * R_FeS2ox
        sediment_tendencies.TCO₂[i, j, 1] = RO₂ + RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄ + R_AOM
        sediment_tendencies.NH₄[i, j, 1] = ratio_NC * (RO₂ +1e-3RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄) + R_DNRA - R_amx - R_NH4ox
        sediment_tendencies.NO₃[i, j, 1] = -2 * RNO₃ - R_DNRA + R_NO2ox
        sediment_tendencies.NO₂[i, j, 1] = 2 * RNO₃ - 1.33 * RNO₂ - R_amx + R_NH4ox - R_NO2ox
        sediment_tendencies.N₂[i, j, 1] = R_amx + 0.66 * RNO₂
        sediment_tendencies.TPO₄[i, j, 1] = ratio_PC * (RO₂ + RNO₃ + RNO₂ + RFe + RSO₄ + RCH₄) + RFe * ratio_FeP - ratio_FeP * R_Fe2ox * fₖ₋ₜₚₒ₄ + ratio_FeP * R_Fe3red # TODO P not conserved???
        sediment_tendencies.FeOHP[i, j, 1] = iron_deposition - 4 * RFe + R_Fe2ox - R_Fe3red
        sediment_tendencies.Feᴵᴵ[i, j, 1] = 4 * RFe - R_Fe2ox + R_FeS2ox - R_FeS2p + R_Fe3red
        sediment_tendencies.FeS₂[i, j, 1] = R_FeS2p - R_FeS2ox
        sediment_tendencies.SO₄[i, j, 1] = -0.5 * RSO₄ + R_DNRA - R_AOM + R_H2Sox + 2 * R_FeS2ox + (R_Fe3red / 8)
        sediment_tendencies.TH₂S[i, j, 1] = 0.5 * RSO₄ - R_DNRA + R_AOM - R_H2Sox  - 2 * R_FeS2p - (R_Fe3red / 8)
        sediment_tendencies.CH₄[i, j, 1] = 0.5 * RCH₄ - R_AOM
        println(sediment_tendencies.FeS₂[i, j, 1])
    end
end

summary(::IronPhosphate{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G + Iron + Phosphate sediment model ($FT)")
show(io::IO, model::IronPhosphate) = print(io, summary(model))
