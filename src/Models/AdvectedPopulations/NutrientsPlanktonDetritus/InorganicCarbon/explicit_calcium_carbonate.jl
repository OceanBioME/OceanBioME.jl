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
using OceanBioME.Models.CarbonChemistryModel: 
    CarbonChemistry, 
    calcium_carbonate_saturation,
    silicate_concentration,
    phosphate_concentration

import Adapt: adapt_structure, adapt
import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    biogeochemical_auxiliary_fields,
    update_biogeochemical_state!

"""
    ExplicitCalciumCarbonate(grid;
                    replicates = 1,
                    calcium_carbonate_dissolution_rate = 0.197 / day,
                    calcium_carbonate_dissolution_exponent = 1,
                    calcium_carbonate_precipitation_rate = 0,
                    calcium_carbonate_precipitation_exponent = 1,
                    calcium_carbonate_sinking_speed = 10 / day,
                    open_bottom = true,
                    carbon_chemistry = CarbonChemistry())

The inorganic carbon component for the `inorganic_carbon` slot of a
[`NutrientsPlanktonDetritus`](@ref) model that tracks calcium carbonate (`CaCO₃`) **explicitly** as its own
sinking tracer, alongside dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`).

Calcifiers form calcium carbonate as they grow (at the plankton's `calcium_carbonate_rain_ratio` per fixed carbon),
removing `DIC` and `2×` that from `Alk`. Each plankton routes its calcium carbonate to three fates through the
`biological_calcium_carbonate_precipitation`/`particulate_calcium_carbonate_production`/`biological_calcium_carbonate_dissolution`
hooks: the calcium carbonate associated with the plankton's **solid** waste enters the explicit, sinking
`CaCO₃` pool; the calcium carbonate associated with its **dissolved and inorganic** waste (e.g. zooplankton gut
dissolution) is returned straight to `DIC`/`Alk`. The `CaCO₃` pool sinks at `calcium_carbonate_sinking_speed`
and dissolves at a saturation-dependent rate

    dissolution = calcium_carbonate_dissolution_rate * max(0, 1 - Ω)^calcium_carbonate_dissolution_exponent * CaCO₃

where `Ω` is the calcium carbonate saturation state (computed each step from `DIC`, `Alk`, `T`, `S` via
`carbon_chemistry`). An optional abiotic precipitation, **off by default**
(`calcium_carbonate_precipitation_rate = 0`), follows a mirror law when the water is supersaturated

    precipitation = calcium_carbonate_precipitation_rate * max(0, Ω - 1)^calcium_carbonate_precipitation_exponent * CaCO₃

Dissolution returns `DIC`/`Alk`; precipitation removes them, so carbon and alkalinity are conserved.
Unlike [`CarbonateSystem`](@ref) (which routes calcium carbonate implicitly through detrital remineralisation),
this component makes `CaCO₃` an explicit pool and so overrides that implicit routing.

`T` and `S` must be present in the model (as for any carbon-chemistry calculation).

Passing `replicates > 1` manifests `replicates` **independent** copies of the carbonate system
(`DIC1`, `Alk1`, `CaCO₃1`, `DIC2`, …); each copy evolves from its own `DIC`/`Alk`/`CaCO₃` and its own
saturation state, which is useful for ensemble or perturbation experiments.

Keyword Arguments
=================
- `grid`: (required) the geometry, needed to build the saturation and sinking-speed fields
- `replicates`: number of independent DIC/Alk/CaCO₃ realisations
- `calcium_carbonate_dissolution_rate`: base dissolution rate `k` (1/s)
- `calcium_carbonate_dissolution_exponent`: exponent `m` of the undersaturation `max(0, 1 - Ω)`
- `calcium_carbonate_precipitation_rate`: abiotic precipitation rate `k_prec` (1/s), `0` disables it
- `calcium_carbonate_precipitation_exponent`: exponent `n` of the supersaturation `max(0, Ω - 1)`
- `calcium_carbonate_sinking_speed`: downward sinking speed of `CaCO₃` (m/s)
- `open_bottom`: whether `CaCO₃` can sink out of the bottom of the domain
- `carbon_chemistry`: the [`CarbonChemistry`](@ref) used to compute the calcium carbonate saturation `Ω`
"""
struct ExplicitCalciumCarbonate{N, FT, CC, SV, SS} <: AbstractInorganicCarbon
          calcium_carbonate_dissolution_rate :: FT
      calcium_carbonate_dissolution_exponent :: FT
        calcium_carbonate_precipitation_rate :: FT
    calcium_carbonate_precipitation_exponent :: FT 
                  carbon_chemistry :: CC
                  sinking_velocity :: SV
                calcium_carbonate_saturation :: SS

    ExplicitCalciumCarbonate{N}(calcium_carbonate_dissolution_rate::FT,
                       calcium_carbonate_dissolution_exponent::FT,
                       calcium_carbonate_precipitation_rate::FT,
                       calcium_carbonate_precipitation_exponent::FT,
                       carbon_chemistry::CC,
                       sinking_velocity::SV,
                       calcium_carbonate_saturation::SS) where {N, FT, CC, SV, SS} =
        new{N, FT, CC, SV, SS}(calcium_carbonate_dissolution_rate,
                               calcium_carbonate_dissolution_exponent,
                               calcium_carbonate_precipitation_rate,
                               calcium_carbonate_precipitation_exponent,
                               carbon_chemistry,
                               sinking_velocity,
                               calcium_carbonate_saturation)
end

function ExplicitCalciumCarbonate(grid::AbstractGrid{FT};
                         replicates = 1,
                         calcium_carbonate_dissolution_rate = convert(FT, 0.197 / day),
                         calcium_carbonate_dissolution_exponent = one(FT),
                         calcium_carbonate_precipitation_rate = zero(FT),
                         calcium_carbonate_precipitation_exponent = one(FT),
                         calcium_carbonate_sinking_speed = convert(FT, 10 / day),
                         open_bottom = true,
                         carbon_chemistry = CarbonChemistry(FT)) where FT

    sinking_velocity = setup_velocity_fields((; CaCO₃ = calcium_carbonate_sinking_speed), grid, open_bottom; three_D = true).CaCO₃

    field_names = saturation_field_names(Val(replicates))
    calcium_carbonate_saturation = NamedTuple{field_names}(ntuple(_ -> CenterField(grid), replicates))

    manifest_explicit_calcium_carbonate_replicates!(replicates)

    return ExplicitCalciumCarbonate{replicates}(convert(FT, calcium_carbonate_dissolution_rate),
                                       convert(FT, calcium_carbonate_dissolution_exponent),
                                       convert(FT, calcium_carbonate_precipitation_rate),
                                       convert(FT, calcium_carbonate_precipitation_exponent),
                                       carbon_chemistry,
                                       sinking_velocity,
                                       calcium_carbonate_saturation)
end

const NPD_EC{FT} = NPD{FT, <:Any, <:Any, <:Any, <:ExplicitCalciumCarbonate}

required_biogeochemical_tracers(::ExplicitCalciumCarbonate{1}) = (:DIC, :Alk, :CaCO₃)
required_biogeochemical_tracers(::ExplicitCalciumCarbonate{N}) where N =
    (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)..., map(n->Symbol(:CaCO₃, n), 1:N)...)

saturation_field_names(::ExplicitCalciumCarbonate{N}) where N = 
    saturation_field_names(Val(N))
saturation_field_names(::Val{1}) = (:Ω, )
saturation_field_names(::Val{N}) where N = ntuple(n -> Symbol(:Ω, n), N)

carbonate_field_names(::ExplicitCalciumCarbonate{1}) = ((:DIC, :Alk), )
carbonate_field_names(::ExplicitCalciumCarbonate{N}) where N = ntuple(n -> (Symbol(:DIC, n), Symbol(:Alk, n)), N)

required_biogeochemical_auxiliary_fields(ic::ExplicitCalciumCarbonate) = saturation_field_names(ic)

@inline function calcium_carbonate_dissolution_flux(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields,
                                          ::Val{calcium_carbonate_name}, ::Val{Ω_name}) where {calcium_carbonate_name, Ω_name}
    Ω     = @inbounds getproperty(auxiliary_fields, Ω_name)[i, j, k]
    CaCO₃ = @inbounds getproperty(fields, calcium_carbonate_name)[i, j, k]

    ic = bgc.inorganic_carbon

    return ic.calcium_carbonate_dissolution_rate * max(zero(Ω), one(Ω) - Ω)^ic.calcium_carbonate_dissolution_exponent * CaCO₃
end

@inline function abiotic_calcium_carbonate_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields,
                                            ::Val{calcium_carbonate_name}, ::Val{Ω_name}) where {calcium_carbonate_name, Ω_name}
    Ω     = @inbounds getproperty(auxiliary_fields, Ω_name)[i, j, k]
    CaCO₃ = @inbounds getproperty(fields, calcium_carbonate_name)[i, j, k]

    ic = bgc.inorganic_carbon

    return ic.calcium_carbonate_precipitation_rate * max(zero(Ω), Ω - one(Ω))^ic.calcium_carbonate_precipitation_exponent * CaCO₃
end

@inline calcium_carbonate_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    biological_calcium_carbonate_precipitation(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) +
    abiotic_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn)

@inline explicit_calcium_carbonate_tendency(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    particulate_calcium_carbonate_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) +
    abiotic_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    calcium_carbonate_dissolution_flux(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn)

@inline net_calcium_carbonate_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields, cn::Val, Ωn::Val) =
    calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    calcium_carbonate_dissolution_flux(i, j, k, grid, bgc, fields, auxiliary_fields, cn, Ωn) -
    biological_calcium_carbonate_dissolution(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

@inline net_calcium_carbonate_production(i, j, k, grid, bgc::NPD_EC, fields, auxiliary_fields) =
    net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, Val(:CaCO₃), Val(:Ω))

@inline (bgc::NPD_EC)(i, j, k, grid, ::Val{:CaCO₃}, clock, fields, auxiliary_fields) =
    explicit_calcium_carbonate_tendency(i, j, k, grid, bgc, fields, auxiliary_fields, Val(:CaCO₃), Val(:Ω))

@inline biogeochemical_drift_velocity(bgc::NPD_EC, ::Val{:CaCO₃}) = bgc.inorganic_carbon.sinking_velocity

biogeochemical_auxiliary_fields(ic::ExplicitCalciumCarbonate) = ic.calcium_carbonate_saturation

function update_biogeochemical_state!(model, ic::ExplicitCalciumCarbonate, npd::NPD)
    for (saturation, carbonate_names) in zip(values(ic.calcium_carbonate_saturation), carbonate_field_names(ic))
        compute_calcium_carbonate_saturation!(ic.carbon_chemistry, saturation, model, carbonate_names)
    end

    return nothing
end

function compute_calcium_carbonate_saturation!(carbon_chemistry, saturation, model, (dic_name, alk_name))
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz, _compute_calcium_carbonate_saturation!,
            carbon_chemistry, saturation, grid, fields(model),
            Val(dic_name), Val(alk_name), defaults.gravitational_acceleration)

    fill_halo_regions!(saturation, model.grid, model.clock, fields(model))

    return nothing
end

@kernel function _compute_calcium_carbonate_saturation!(carbon_chemistry, saturation, grid, model_fields,
                                              ::Val{dic_name}, ::Val{alk_name}, g) where {dic_name, alk_name}
    i, j, k = @index(Global, NTuple)

    T   = @inbounds model_fields.T[i, j, k]
    S   = @inbounds model_fields.S[i, j, k]
    DIC = @inbounds getproperty(model_fields, dic_name)[i, j, k]
    Alk = @inbounds getproperty(model_fields, alk_name)[i, j, k]

    silicate  = silicate_concentration(grid, i, j, k, model_fields)
    phosphate = phosphate_concentration(grid, i, j, k, model_fields)

    z = znode(i, j, k, grid, Center(), Center(), Center())

    P = abs(z) * g * 1026 / 100000

    @inbounds saturation[i, j, k] = calcium_carbonate_saturation(carbon_chemistry; DIC, T, S, Alk, P, phosphate, silicate)
end

const _manifested_explicit_calcium_carbonate = Set{Int}()

function manifest_explicit_calcium_carbonate_replicates!(N)
    (N > 1 && !(N in _manifested_explicit_calcium_carbonate)) || return nothing
    push!(_manifested_explicit_calcium_carbonate, N)

    for n in 1:N
        DIC_name   = Symbol(:DIC, n)
        Alk_name   = Symbol(:Alk, n)
        CaCO₃_name = Symbol(:CaCO₃, n)
        Ω_name     = Symbol(:Ω, n)
        @eval begin
            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(CaCO₃_name))}, clock, fields, auxiliary_fields) =
                explicit_calcium_carbonate_tendency(i, j, k, grid, bgc, fields, auxiliary_fields,
                                          Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) -
                net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields,
                                       Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline (bgc::NPD_EC)(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields) -
                2 * net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields,
                                           Val($(QuoteNode(CaCO₃_name))), Val($(QuoteNode(Ω_name))))

            @inline biogeochemical_drift_velocity(bgc::NPD_EC, ::Val{$(QuoteNode(CaCO₃_name))}) =
                bgc.inorganic_carbon.sinking_velocity
        end
    end

    return nothing
end

Adapt.adapt_structure(to, ic::ExplicitCalciumCarbonate{N}) where N =
    ExplicitCalciumCarbonate{N}(ic.calcium_carbonate_dissolution_rate,
                       ic.calcium_carbonate_dissolution_exponent,
                       ic.calcium_carbonate_precipitation_rate,
                       ic.calcium_carbonate_precipitation_exponent,
                       adapt(to, ic.carbon_chemistry),
                       adapt(to, ic.sinking_velocity),
                       adapt(to, ic.calcium_carbonate_saturation))

Base.summary(ic::ExplicitCalciumCarbonate{1}) =
    string("ExplicitCalciumCarbonate $(required_biogeochemical_tracers(ic))")

Base.summary(ic::ExplicitCalciumCarbonate{N}) where N =
    string("ExplicitCalciumCarbonate{realisations = $N} $(required_biogeochemical_tracers(ic))")

function Base.show(io::IO, ic::ExplicitCalciumCarbonate{N}) where N
    msg = "ExplicitCalciumCarbonate $(required_biogeochemical_tracers(ic))"

    if N > 1
        msg *= "\n└── Realisations: $N"
    end

    print(io, msg)

    return nothing
end
