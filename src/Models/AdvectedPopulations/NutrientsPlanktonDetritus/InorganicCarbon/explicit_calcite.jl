using Adapt
using Oceananigans: defaults
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: CenterField
using Oceananigans.Grids: AbstractGrid, znode, Center
using Oceananigans.Models: fields
using Oceananigans.Units: day
using Oceananigans.Utils: launch!

using KernelAbstractions: @kernel, @index

using OceanBioME: setup_velocity_fields
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry, calcite_saturation

import Adapt: adapt_structure, adapt
import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    biogeochemical_auxiliary_fields,
    update_biogeochemical_state!

"""
    ExplicitCalcite(grid;
                    replicates = 1,
                    calcite_dissolution_rate = 0.197 / day,
                    calcite_dissolution_exponent = 1,
                    calcite_precipitation_rate = 0,
                    calcite_precipitation_exponent = 1,
                    calcite_sinking_speed = 10 / day,
                    open_bottom = true,
                    carbon_chemistry = CarbonChemistry())

The inorganic carbon component for the `inorganic_carbon` slot of a
[`NutrientsPlanktonDetritus`](@ref) model that tracks calcite (`CaCO₃`) **explicitly** as its own
sinking tracer, alongside dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`).

Calcifiers form calcite as they grow (at the plankton's `calcite_rain_ratio` per fixed carbon),
removing `DIC` and `2×` that from `Alk`. Each plankton routes its calcite to three fates through the
`biological_calcite_precipitation`/`particulate_calcite_production`/`biological_calcite_dissolution`
hooks: the calcite associated with the plankton's **solid** waste enters the explicit, sinking
`CaCO₃` pool; the calcite associated with its **dissolved and inorganic** waste (e.g. zooplankton gut
dissolution) is returned straight to `DIC`/`Alk`. The `CaCO₃` pool sinks at `calcite_sinking_speed`
and dissolves at a saturation-dependent rate

    dissolution = calcite_dissolution_rate * max(0, 1 - Ω)^calcite_dissolution_exponent * CaCO₃

where `Ω` is the calcite saturation state (computed each step from `DIC`, `Alk`, `T`, `S` via
`carbon_chemistry`). An optional abiotic precipitation, **off by default**
(`calcite_precipitation_rate = 0`), follows a mirror law when the water is supersaturated

    precipitation = calcite_precipitation_rate * max(0, Ω - 1)^calcite_precipitation_exponent * CaCO₃

Dissolution returns `DIC`/`Alk`; precipitation removes them, so carbon and alkalinity are conserved.
Unlike [`CarbonateSystem`](@ref) (which routes calcite implicitly through detrital remineralisation),
this component makes `CaCO₃` an explicit pool and so overrides that implicit routing.

`T` and `S` must be present in the model (as for any carbon-chemistry calculation).

Passing `replicates > 1` manifests `replicates` **independent** copies of the carbonate system
(`DIC1`, `Alk1`, `CaCO₃1`, `DIC2`, …); each copy evolves from its own `DIC`/`Alk`/`CaCO₃` and its own
saturation state, which is useful for ensemble or perturbation experiments.

Keyword Arguments
=================
- `grid`: (required) the geometry, needed to build the saturation and sinking-speed fields
- `replicates`: number of independent DIC/Alk/CaCO₃ realisations
- `calcite_dissolution_rate`: base dissolution rate `k` (1/s)
- `calcite_dissolution_exponent`: exponent `m` of the undersaturation `max(0, 1 - Ω)`
- `calcite_precipitation_rate`: abiotic precipitation rate `k_prec` (1/s), `0` disables it
- `calcite_precipitation_exponent`: exponent `n` of the supersaturation `max(0, Ω - 1)`
- `calcite_sinking_speed`: downward sinking speed of `CaCO₃` (m/s)
- `open_bottom`: whether `CaCO₃` can sink out of the bottom of the domain
- `carbon_chemistry`: the [`CarbonChemistry`](@ref) used to compute the calcite saturation `Ω`
"""
struct ExplicitCalcite{N, SN, FT, CC, SV, ST} <: AbstractInorganicCarbon
          calcite_dissolution_rate :: FT
      calcite_dissolution_exponent :: FT
        calcite_precipitation_rate :: FT
    calcite_precipitation_exponent :: FT 
                  carbon_chemistry :: CC
                  sinking_velocity :: SV
                calcite_saturation :: NamedTuple{SN, ST} 

    ExplicitCalcite{N}(calcite_dissolution_rate::FT,
                       calcite_dissolution_exponent::FT,
                       calcite_precipitation_rate::FT,
                       calcite_precipitation_exponent::FT,
                       carbon_chemistry::CC,
                       sinking_velocity::SV,
                       calcite_saturation::NamedTuple{SN, ST}) where {N, SN, FT, CC, SV, ST} =
        new{N, SN, FT, CC, SV, ST}(calcite_dissolution_rate,
                                   calcite_dissolution_exponent,
                                   calcite_precipitation_rate,
                                   calcite_precipitation_exponent,
                                   carbon_chemistry,
                                   sinking_velocity,
                                   calcite_saturation)
end

function ExplicitCalcite(grid::AbstractGrid{FT};
                         replicates = 1,
                         calcite_dissolution_rate = convert(FT, 0.197 / day),
                         calcite_dissolution_exponent = one(FT),
                         calcite_precipitation_rate = zero(FT),
                         calcite_precipitation_exponent = one(FT),
                         calcite_sinking_speed = convert(FT, 10 / day),
                         open_bottom = true,
                         carbon_chemistry = CarbonChemistry(FT)) where FT

    sinking_velocity = setup_velocity_fields((; CaCO₃ = calcite_sinking_speed), grid, open_bottom; three_D = true).CaCO₃

    field_names = saturation_field_names(Val(replicates))
    calcite_saturation = NamedTuple{field_names}(ntuple(_ -> CenterField(grid), replicates))

    manifest_explicit_calcite_replicates!(replicates)

    return ExplicitCalcite{replicates}(convert(FT, calcite_dissolution_rate),
                                       convert(FT, calcite_dissolution_exponent),
                                       convert(FT, calcite_precipitation_rate),
                                       convert(FT, calcite_precipitation_exponent),
                                       carbon_chemistry,
                                       sinking_velocity,
                                       calcite_saturation)
end

const NPD_EC{FT} = NPD{FT, <:Any, <:Any, <:Any, <:ExplicitCalcite}

required_biogeochemical_tracers(::ExplicitCalcite{1}) = (:DIC, :Alk, :CaCO₃)
required_biogeochemical_tracers(::ExplicitCalcite{N}) where N =
    (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)..., map(n->Symbol(:CaCO₃, n), 1:N)...)

# per-realisation saturation-field names, and the (DIC, Alk) pairs each is computed from
saturation_field_names(::Val{1}) = (:Ω, )
saturation_field_names(::Val{N}) where N = ntuple(n -> Symbol(:Ω, n), N)

carbonate_field_names(::ExplicitCalcite{1}) = ((:DIC, :Alk), )
carbonate_field_names(::ExplicitCalcite{N}) where N = ntuple(n -> (Symbol(:DIC, n), Symbol(:Alk, n)), N)

required_biogeochemical_auxiliary_fields(::ExplicitCalcite{N, FN}) where {N, FN} = FN

@inline function calcite_dissolution_flux(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields,
                                          ::Val{calcite_name}, ::Val{Ω_name}) where {calcite_name, Ω_name}
    Ω     = @inbounds getproperty(auxiliary_fields, Ω_name)[i, j, k]
    CaCO₃ = @inbounds getproperty(fields, calcite_name)[i, j, k]

    ic = bgc.inorganic_carbon

    return ic.calcite_dissolution_rate * max(zero(Ω), one(Ω) - Ω)^ic.calcite_dissolution_exponent * CaCO₃
end

@inline function abiotic_calcite_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields,
                                            ::Val{calcite_name}, ::Val{Ω_name}) where {calcite_name, Ω_name}
    Ω     = @inbounds getproperty(auxiliary_fields, Ω_name)[i, j, k]
    CaCO₃ = @inbounds getproperty(fields, calcite_name)[i, j, k]

    ic = bgc.inorganic_carbon

    return ic.calcite_precipitation_rate * max(zero(Ω), Ω - one(Ω))^ic.calcite_precipitation_exponent * CaCO₃
end

@inline calcite_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    biological_calcite_precipitation(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) +
    abiotic_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn)

@inline explicit_calcite_tendency(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    particulate_calcite_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) +
    abiotic_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    calcite_dissolution_flux(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn)

@inline net_calcite_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    calcite_dissolution_flux(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    biological_calcite_dissolution(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

@inline net_calcite_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields) =
    net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields, Val(:CaCO₃), Val(:Ω))

@inline (bgc::NPD_EC)(i, j, k, grid, ::Val{:CaCO₃}, clock, fields, auxiliary_fields) =
    explicit_calcite_tendency(i, j, k, grid, bgc, fields, auxiliary_fields, Val(:CaCO₃), Val(:Ω))

#####
##### sinking
#####

@inline biogeochemical_drift_velocity(bgc::NPD_EC, ::Val{:CaCO₃}) = bgc.inorganic_carbon.sinking_velocity

#####
##### calcite saturation Ω — precomputed once per step (avoids repeating the pH Newton solve per RK
##### substep). Requires T and S to be present, as for any carbon-chemistry calculation.
#####

biogeochemical_auxiliary_fields(ic::ExplicitCalcite) = ic.calcite_saturation

function update_biogeochemical_state!(model, ic::ExplicitCalcite, npd::NPD)
    for (saturation, carbonate_names) in zip(ic.calcite_saturation, carbonate_field_names(ic))
        compute_calcite_saturation!(ic.carbon_chemistry, saturation, model, carbonate_names)
    end

    return nothing
end

# alias to reference the CarbonChemistry function inside code that shadows `calcite_saturation`
const calcite_saturation_state = calcite_saturation

function compute_calcite_saturation!(carbon_chemistry, saturation, model, (dic_name, alk_name))
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz, _compute_calcite_saturation!,
            carbon_chemistry, saturation, grid, fields(model),
            Val(dic_name), Val(alk_name), defaults.gravitational_acceleration)

    fill_halo_regions!(saturation, model.grid, model.clock, fields(model))

    return nothing
end

@kernel function _compute_calcite_saturation!(carbon_chemistry, saturation, grid, model_fields,
                                              ::Val{dic_name}, ::Val{alk_name}, g) where {dic_name, alk_name}
    i, j, k = @index(Global, NTuple)

    T   = @inbounds model_fields.T[i, j, k]
    S   = @inbounds model_fields.S[i, j, k]
    DIC = @inbounds getproperty(model_fields, dic_name)[i, j, k]
    Alk = @inbounds getproperty(model_fields, alk_name)[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    # rough hydrostatic pressure (bar); matches the PISCES calcite-saturation approximation
    P = abs(z) * g * 1026 / 100000

    @inbounds saturation[i, j, k] = calcite_saturation_state(carbon_chemistry; DIC, T, S, Alk, P)
end

const _manifested_explicit_calcite = Set{Int}()

function manifest_explicit_calcite_replicates!(N)
    (N > 1 && !(N in _manifested_explicit_calcite)) || return nothing
    push!(_manifested_explicit_calcite, N)

    for n in 1:N
        DIC_name   = Symbol(:DIC, n)
        Alk_name   = Symbol(:Alk, n)
        CaCO₃_name = Symbol(:CaCO₃, n)
        Ω_name     = Symbol(:Ω, n)
        @eval begin
            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(CaCO₃_name))}, clock, fields, auxiliary_fields) =
                explicit_calcite_tendency(i, j, k, grid, bgc, fields, auxiliary_fields,
                                          Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) -
                net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields,
                                       Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields) -
                2 * net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields,
                                           Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline biogeochemical_drift_velocity(bgc::NPD_EC, ::Val{$(QuoteNode(CaCO₃_name))}) =
                bgc.inorganic_carbon.sinking_velocity
        end
    end

    return nothing
end

#####
##### adapt / show
#####

Adapt.adapt_structure(to, ic::ExplicitCalcite{N}) where N =
    ExplicitCalcite{N}(ic.calcite_dissolution_rate,
                       ic.calcite_dissolution_exponent,
                       ic.calcite_precipitation_rate,
                       ic.calcite_precipitation_exponent,
                       adapt(to, ic.carbon_chemistry),
                       adapt(to, ic.sinking_velocity),
                       adapt(to, ic.calcite_saturation))

Base.summary(ic::ExplicitCalcite{1}) =
    string("ExplicitCalcite $(required_biogeochemical_tracers(ic))")

Base.summary(ic::ExplicitCalcite{N}) where N =
    string("ExplicitCalcite{realisations = $N} $(required_biogeochemical_tracers(ic))")

function Base.show(io::IO, ic::ExplicitCalcite{N}) where N
    msg = "ExplicitCalcite $(required_biogeochemical_tracers(ic))"

    if N > 1
        msg *= "\n└── Realisations: $N"
    end

    print(io, msg)

    return nothing
end
