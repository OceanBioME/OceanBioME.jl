using Oceananigans.Units

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Fields: Center
using Oceananigans.Grids: znode

import Base: show, summary

"""
    struct SimpleMultiG

Hold the parameters and fields for a simple "multi G" single-layer sediment model.
Based on the Level 3 model described by [Soetaert2000](@citet).
"""
struct SimpleMultiG{SR, FT, P1, P2, P3, P4, SN, SC} <: AbstractContinuousFormSedimentBiogeochemistry
            sinking_redfield :: SR

             fast_decay_rate :: FT
             slow_decay_rate :: FT

               fast_redfield :: FT
               slow_redfield :: FT

               fast_fraction :: FT
               slow_fraction :: FT
          refactory_fraction :: FT

          sedimentation_rate :: FT
      anoxia_half_saturation :: FT

    nitrate_oxidation_params :: P1
      denitrification_params :: P2
               anoxic_params :: P3
            solid_dep_params :: P4

            sinking_nitrogen :: SN
              sinking_carbon :: SC
end

@inline required_sediment_fields(::SimpleMultiG{Nothing}) = (:Cs, :Cf, :Cr, :Ns, :Nf, :Nr)
@inline required_sediment_fields(::SimpleMultiG) = (:Ns, :Nf, :Nr)
@inline required_tracers(::SimpleMultiG) = (:NO₃, :NH₄, :O₂)
@inline sinking_fluxs(s::SimpleMultiG{Nothing}) = (s.sinking_nitrogen..., s.sinking_carbon...)
@inline sinking_fluxs(s::SimpleMultiG) = s.sinking_nitrogen
@inline coupled_tracers(::SimpleMultiG) =  (:NO₃, :NH₄, :O₂)
@inline coupled_tracers(::SimpleMultiG{Nothing}) =  (:NO₃, :NH₄, :O₂, :DIC)

"""
    SimpleMultiGSediment(grid;
                          ...)

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

julia> sediment_model = SimpleMultiG(grid)
...
```
"""
SimpleMultiGSediment(grid;
                     fast_decay_rate = 2/day,
                     slow_decay_rate = 0.2/day,
                     fast_redfield = 0.1509,
                     slow_redfield = 0.13,
                     fast_fraction = 0.74,
                     slow_fraction = 0.26,
                     refactory_fraction = 0.1,
                     sedimentation_rate = 982 * abs(znode(1, 1, 1, grid, Center(), Center(), Center())) ^ (-1.548), # cm/year, incorrect for D < 100m
                     anoxia_half_saturation = 1.0, # mmol/m³ (arbitarily low)
                     nitrate_oxidation_params = on_architecture(architecture(grid), (- 1.9785, 0.2261, -0.0615, -0.0289, - 0.36109, - 0.0232)),
                     denitrification_params = on_architecture(architecture(grid), (- 3.0790, 1.7509, 0.0593, - 0.1923, 0.0604, 0.0662)),
                     anoxic_params = on_architecture(architecture(grid), (- 3.9476, 2.6269, - 0.2426, -1.3349, 0.1826, - 0.0143)),
                     solid_dep_params = on_architecture(architecture(grid), (0.233, 0.336, 982.0, - 1.548)),
                     sinking_nitrogen = (:sPOM, :bPOM),
                     sinking_carbon = nothing,
                     sinking_redfield = ifelse(isnothing(sinking_carbon), convert(eltype(grid), 6.56), nothing),
                     kwargs...) =
    BiogeochemicalSediment(grid, 
                           SimpleMultiG(sinking_redfield, 
                                        fast_decay_rate, slow_decay_rate,
                                        fast_redfield,   slow_redfield,
                                        fast_fraction,   slow_fraction, refactory_fraction,
                                        sedimentation_rate, anoxia_half_saturation,
                                        nitrate_oxidation_params, denitrification_params,
                                        anoxic_params, solid_dep_params,
                                        tracernames(sinking_nitrogen),
                                        tracernames(sinking_carbon));
                           kwargs...)

adapt_structure(to, sediment::SimpleMultiG) = 
    SimpleMultiG(adapt(to, sediment.sinking_redfield),
                 adapt(to, sediment.fast_decay_rate),
                 adapt(to, sediment.slow_decay_rate),
                 adapt(to, sediment.fast_redfield),
                 adapt(to, sediment.slow_redfield),
                 adapt(to, sediment.fast_fraction),
                 adapt(to, sediment.slow_fraction),
                 adapt(to, sediment.refactory_fraction),
                 adapt(to, sediment.sedimentation_rate),
                 adapt(to, sediment.anoxia_half_saturation),
                 adapt(to, sediment.nitrate_oxidation_params),
                 adapt(to, sediment.denitrification_params),
                 adapt(to, sediment.anoxic_params),
                 adapt(to, sediment.solid_dep_params),
                 adapt(to, sediment.sinking_nitrogen),
                 adapt(to, sediment.sinking_carbon))

@inline sinking_nitrogen_carbon(PON, POC) = PON, POC
@inline sinking_nitrogen_carbon(sPON, bPON, sPOC, bPOC) = sPON+bPON, sPOC+bPOC

@inline function sinking_nitrogen_carbon(args...)
    n = Int(length(args)/2)

    N = @inbounds sum(args[1:n])
    C = @inbounds sum(args[n+1:end])

    return N, C
end

### nitrogen only

@inline function (s::SimpleMultiG)(::Val{:Ns}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_nitrogen...)
    λ = s.slow_decay_rate
    fr = s.refactory_fraction
    fs = s.slow_fraction

    flux = sum(sinking_nitrogen)

    return (1 - fr) * fs * flux - λ * Ns
end

@inline function (s::SimpleMultiG)(::Val{:Nf}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_nitrogen...)
    λ = s.fast_decay_rate
    fr = s.refactory_fraction
    ff = s.fast_fraction

    flux = sum(sinking_nitrogen)

    return (1 - fr) * ff * flux - λ * Nf
end

@inline function (s::SimpleMultiG)(::Val{:Nr}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_nitrogen...)
    fr = s.refactory_fraction

    flux = sum(sinking_nitrogen)

    return fr * flux
end

# tracer tendencies

@inline function (s::SimpleMultiG)(::Val{:NH₄}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_tracers...)
    R  = s.sinking_redfield
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    Nr = λs*Ns + λf*Nf
    Cr = Nr * R  

    k = reactivity(s, Ns * R, Nf * R) 

    pₙ  = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)

    return (1 - pₙ) * Nr + 0.8 * pₙ′ * Cr # I think the implication is that this `0.8 * pₙ′ * Cr` term should be going to N₂, but we want to conserve
end

@inline function (s::SimpleMultiG)(::Val{:NO₃}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_tracers...)
    R  = s.sinking_redfield
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    Nr = λs*Ns + λf*Nf
    Cr = Nr * R  

    k = reactivity(s, Ns * R, Nf * R) 

    pₙ  = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)

    # when pₙ > 1 there is some instant conversion of NH₄ to NO₃ which is kind of strange
    return pₙ * Nr - 0.8 * pₙ′ * Cr
end

@inline function (s::SimpleMultiG)(::Val{:O₂}, x, y, t, Ns, Nf, Nr, NO₃, NH₄, O₂, sinking_tracers...)
    R  = s.sinking_redfield
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    kO₂ = s.anoxia_half_saturation

    Nr = λs*Ns + λf*Nf
    Cr = Nr * R  

    k = reactivity(s, Ns * R, Nf * R) 

    pₙ  = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)
    pₐ  = anoxic_remineralisation_fraction(s, Nr, Cr, k, NO₃, O₂)

    pₛ  = solid_deposition_fraction(s)

    return - (1 - pₐ * pₛ - pₙ′) * O₂ / (kO₂ + O₂) * Cr - 2pₙ * Nr # here there is an O:C of 1 and O:N of 2
end

### nitrogen and carbon

@inline function (s::SimpleMultiG{Nothing})(::Val{:Ns}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λ = s.slow_decay_rate
    fr = s.refactory_fraction
    fs = s.slow_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return (1 - fr) * fs * fN - λ * Ns
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:Nf}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λ = s.fast_decay_rate
    fr = s.refactory_fraction
    ff = s.fast_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return (1 - fr) * ff * fN - λ * Nf
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:Nr}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    fr = s.refactory_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return fr * fN
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:Cs}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λ = s.slow_decay_rate
    fr = s.refactory_fraction
    fs = s.slow_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return (1 - fr) * fs * fC - λ * Cs
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:Cf}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λ = s.fast_decay_rate
    fr = s.refactory_fraction
    ff = s.fast_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return (1 - fr) * ff * fC - λ * Cf
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:Cr}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    fr = s.refactory_fraction

    fN, fC = sinking_nitrogen_carbon(sinking_tracers...)

    return fr * fC
end


# tracer tendencies

@inline function (s::SimpleMultiG{Nothing})(::Val{:NH₄}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    Nr = λs*Ns + λf*Nf
    Cr = λs*Cs + λf*Cf

    k = reactivity(s, Cs, Cf) 

    pₙ = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)

    return (1 - pₙ) * Nr + 0.8 * pₙ′ * Cr # I think the implication is that this `0.8 * pₙ′ * Cr` term should be going to N₂, but we want to conserve
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:NO₃}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    Nr = λs*Ns + λf*Nf
    Cr = λs*Cs + λf*Cf

    k = reactivity(s, Cs, Cf) 

    pₙ  = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)

    # when pₙ > 1 there is some instant conversion of NH₄ to NO₃ which is kind of strange
    return pₙ * Nr - 0.8 * pₙ′ * Cr
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:O₂}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    kO₂ = s.anoxia_half_saturation

    Nr = λs*Ns + λf*Nf
    Cr = λs*Cs + λf*Cf

    k = reactivity(s, Cs, Cf) 

    pₙ  = ammonia_oxidation_fraction(s, Nr, Cr, k, NH₄, O₂)
    pₙ′ = denitrifcation_fraction(s, Nr, Cr, k, NO₃, O₂)
    pₐ  = anoxic_remineralisation_fraction(s, Nr, Cr, k, NO₃, O₂)

    pₛ  = solid_deposition_fraction(s)

    return - (1 - pₐ * pₛ - pₙ′) * O₂ / (kO₂ + O₂) * Cr - 2pₙ * Nr # here there is an O:C of 1 and O:N of 2
end

@inline function (s::SimpleMultiG{Nothing})(::Val{:DIC}, x, y, t, Ns, Nf, Nr, Cs, Cf, Cr, NO₃, NH₄, O₂, sinking_tracers...)
    λs = s.slow_decay_rate
    λf = s.fast_decay_rate

    Cr = λs*Cs + λf*Cf

    return Cr
end

### reactivity constants

@inline function reactivity(sediment, Cs, Cf)
    λs = sediment.slow_decay_rate
    λf = sediment.fast_decay_rate

    Cr = λs*Cs + λf*Cf

    return Cr / (Cs + Cf + eps(0.0))
end

@inline function ammonia_oxidation_fraction(sediment, Nr, Cr, k, NH₄, O₂)
    A, B, C, D, E, F = sediment.nitrate_oxidation_params

    kO₂ = sediment.anoxia_half_saturation

    ln_pNr = (A + B * log(Cr * day) * log(O₂) 
                + C * log(Cr * day)^2 
                + D * log(k * day) * log(NH₄) 
                + E * log(Cr * day) 
                + F * log(Cr * day) * log(NH₄))

    p = exp(ln_pNr) / (Nr * day) * O₂ / (kO₂ + O₂)

    return ifelse(isfinite(p), p, zero(Nr))
end

@inline function denitrifcation_fraction(sediment, Nr, Cr, k, NO₃, O₂)
    A, B, C, D, E, F = sediment.denitrification_params

    kO₂ = sediment.anoxia_half_saturation

    ln_pCr = (A + B * log(Cr * day) 
                + C * log(NO₃)^2 
                + D * log(Cr * day)^2 
                + E * log(k * day)^2
                + F * log(O₂) * log(k))

    p = exp(ln_pCr) / (Cr * day) * O₂ / (kO₂ + O₂)

    return ifelse(isfinite(p), p, zero(Nr))
end

@inline function anoxic_remineralisation_fraction(sediment, Nr, Cr, k, NO₃, O₂)
    A, B, C, D, E, F = sediment.anoxic_params

    ln_pCr = (A + B * log(Cr * day) 
                + C * log(Cr * day)^2 
                + D * log(k * day) 
                + E * log(O₂) * log(k)
                + F * log(NO₃)^2)

    p = exp(ln_pCr) / (Cr * day)

    return ifelse(isfinite(p), p, zero(Nr))
end

@inline solid_deposition_fraction(s) = 0.223 * s.sedimentation_rate ^ 0.336

summary(::SimpleMultiG{FT, P1, P2, P3, P4, F, TE}) where {FT, P1, P2, P3, P4, F, TE} = string("Single-layer multi-G sediment model ($FT)")
show(io::IO, model::SimpleMultiG) = print(io, summary(model))
