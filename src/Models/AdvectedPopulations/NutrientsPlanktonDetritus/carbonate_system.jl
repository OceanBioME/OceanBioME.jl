using Adapt

using Oceananigans: defaults, Units
using Oceananigans.Architectures: architecture
using Oceananigans.Fields: CenterField, Center
using Oceananigans.Grids: znode, AbstractGrid, Flat
using Oceananigans.Utils: launch!

using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry, calcite_saturation

using KernelAbstractions: @kernel, @index

import Adapt: adapt_structure

"""
    CarbonateSystem{N, CAL}

`CarbonateSystem` defines the carbonate system for the `NutrientsPlanktonDetritus` biogeochemical
model and evolves dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`).
`DIC` is the total dissolved inorganic carbon, typically considered to be dissolved 
carbon dioxide, carbonic acid, bicarbonate, and carbonate. `Alk`alinity is 
the buffering capacity of the water which is approximately the amount of bases
in solution that can be neutralised by the addition of acid so is approximately
the sum of non-conserved ions such as carbonate, bicarbonate, hydroxide ions,
etc.

With this module activated the carbon budget can be closed as the flux of carbon
into the plankton and detritus comes from the DIC. Additionally, the flux in/out 
of the water from the air can be computed via the `CarbonChemistry` module.

`DIC` is produced by remineralisation of organic matter, for example the breakdown
of detritus or excretion from zooplankton waste. It is also produced by the 
dissolution of calcite in the zooplankton gut. It is removed by photosynthesis
in phytoplankton.

`Alk`alinity is modified by nitrogen cycling and calcification. Nitrate uptake by phytoplankton increases alkalinity, while calcification decreases alkalinity through the formation of CaCO₃.

`DIC` and `Alk` concentration is only one way coupled with the rest of the 
biogeochemistry and *does not* effect any other groups (e.g. acidification does
not affect phytoplankton growth). To capture this effect a different `plankton` 
could be defined.

Multiple (N) instances of the carbonate system can evolve in parallel and the tracers
will be named `DIC1` ... `DICN` and `Alk1` ... `AlkN`. This is setup by passing an 
integer argument when the model is constructed like `CarbonateSystem(2)`.

## Explicit calcite

Pass a `grid` and set `explicit_calcite = true` to evolve a `CaCO₃` tracer with kinetic
precipitation and dissolution. This requires temperature (`T`, °C) and salinity (`S`, PSU)
to be model tracers so that the calcite saturation state (Ω) can be computed each timestep.
Include them when constructing the ocean model, e.g.
`NonhydrostaticModel(grid; tracers = (:T, :S), biogeochemistry = NutrientsPlanktonDetritus(grid; carbonate_system = CarbonateSystem(grid; explicit_calcite=true)))`.

You may get method overwrite warnings if you repeatedly define carbonate systems
with N>1, but this shouldn't cause problems.
"""
struct CarbonateSystem{N, CAL, CSS, CC, FT, SV}
      calcite_saturation_state :: CSS
              carbon_chemistry :: CC
    gravitational_acceleration :: FT
     dissolution_rate_constant :: FT 
             dissolution_order :: FT
   precipitation_rate_constant :: FT
           precipitation_order :: FT
      calcite_sinking_velocity :: SV

    CarbonateSystem{N, CAL}(calcite_saturation_state::CSS,
                            carbon_chemistry::CC,
                            gravitational_acceleration::FT,
                            dissolution_rate_constant::FT,
                            dissolution_order::FT,
                            precipitation_rate_constant::FT,
                            precipitation_order::FT,
                            calcite_sinking_velocity::SV) where {N, CAL, CSS, CC, FT, SV} =
        new{N, CAL, CSS, CC, FT, SV}(calcite_saturation_state,
                                     carbon_chemistry,
                                     gravitational_acceleration,
                                     dissolution_rate_constant,
                                     dissolution_order,
                                     precipitation_rate_constant,
                                     precipitation_order,
                                     calcite_sinking_velocity)
end

function CarbonateSystem(replicates::Int = 1)
    manifest_carbonate_replicates!(replicates)

    return CarbonateSystem{replicates, false}(nothing, nothing, nothing, nothing,
                                              nothing, nothing, nothing, nothing)
end

"""
    CarbonateSystem(grid, replicates = 1; explicit_calcite = false, ...)

Construct a [`CarbonateSystem`](@ref) on `grid`, optionally with an explicit `CaCO₃` tracer.

When `explicit_calcite = true`, the ocean model must include `T` and `S` tracers (°C and PSU) for calcite saturation calculations. These are included in `required_biogeochemical_tracers` alongside `DIC`, `Alk`, and `CaCO₃`.
"""
function CarbonateSystem(grid::AbstractGrid, replicates = 1; 
                         explicit_calcite = false,
                         carbon_chemistry::CC = CarbonChemistry(),
                         gravitational_acceleration = defaults.gravitational_acceleration,
                         dissolution_rate_constant = 1e-4 / day,
                         dissolution_order = 1.0,
                         precipitation_rate_constant = 1e-5 / day,
                         precipitation_order = 2.0,
                         sinking_speed = 10 / day,
                         open_bottom = true) where CC
    
    if !explicit_calcite
        return CarbonateSystem(replicates)
    end

    manifest_carbonate_replicates!(replicates)
    
    calcite_saturation_state = tuple(map(n->CenterField(grid), 1:replicates)...)
    sinking_velocities = setup_velocity_fields((; C = sinking_speed), grid, open_bottom; three_D = true).C

    return CarbonateSystem{replicates, explicit_calcite}(calcite_saturation_state, 
                                                         carbon_chemistry, 
                                                         gravitational_acceleration,
                                                         dissolution_rate_constant,
                                                         dissolution_order,
                                                         precipitation_rate_constant,
                                                         precipitation_order,
                                                         sinking_velocities)
end

Adapt.adapt_structure(to, cs::CarbonateSystem{N, CAL}) where {N, CAL} = 
    CarbonateSystem{N, CAL}(adapt(to, cs.calcite_saturation_state),
                            adapt(to, cs.carbon_chemistry),
                            adapt(to, cs.gravitational_acceleration),
                            adapt(to, cs.dissolution_rate_constant),
                            adapt(to, cs.dissolution_order),
                            adapt(to, cs.precipitation_rate_constant),
                            adapt(to, cs.precipitation_order),
                            adapt(to, cs.calcite_sinking_velocity))

function update_biogeochemical_state!(model, bgc::CarbonateSystem{N, true}) where N
    carbon_chemistry = bgc.carbon_chemistry
    grid = model.grid

    # T and S must be model tracers when explicit_calcite = true; see required_biogeochemical_tracers
    T = model.tracers.T
    S = model.tracers.S

    arch = architecture(grid)

    for n in 1:N
        DIC = model.tracers[Symbol(:DIC, ifelse(N==1, "", n))]
        Alk = model.tracers[Symbol(:Alk, ifelse(N==1, "", n))]

        Ω   = @inbounds bgc.calcite_saturation_state[n]

        launch!(arch, model.grid, :xyz, 
                compute_calcite_saturation!, 
                Ω, DIC, Alk, T, S, grid, carbon_chemistry, bgc)
    end
end

@inline function pressure(i, j, k, grid::AbstractGrid{FT}, g, ρ) where FT
    z = znode(i, j, k, grid, Center(), Center(), Center())

    return abs(z) * convert(FT, ρ / 1e5) * g
end

@inline pressure(i, j, k, grid::AbstractGrid{<:Any, Flat, Flat, Flat}, g, ρ) = zero(grid)

@kernel function compute_calcite_saturation!(Ω, DIC, Alk, T, S, grid::AbstractGrid{FT}, carbon_chemistry, bgc) where FT
    i, j, k = @index(Global, NTuple)
    
    g = bgc.gravitational_acceleration

    P = pressure(i, j, k, grid, g, 1026)

    @inbounds begin
        Ω[i, j, k] = calcite_saturation(carbon_chemistry;
                                        T = T[i, j, k],
                                        S = S[i, j, k],
                                        DIC = DIC[i, j, k],
                                        Alk = Alk[i, j, k],
                                        P)
    end
end

required_biogeochemical_tracers(::CarbonateSystem{1, false}) = (:DIC, :Alk)
required_biogeochemical_tracers(::CarbonateSystem{1, true}) = (:DIC, :Alk, :CaCO₃, :T, :S)

required_biogeochemical_tracers(::CarbonateSystem{N, false}) where N = (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)
required_biogeochemical_tracers(::CarbonateSystem{N, true}) where N = (map(n->Symbol(:DIC, n), 1:N)..., 
                                                                       map(n->Symbol(:Alk, n), 1:N)...,
                                                                       map(n->Symbol(:CaCO₃, n), 1:N)...,
                                                                       :T, :S)

@inline biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, 
                                      ::Val{:CaCO₃}) = 
    bgc.carbonate_system.calcite_sinking_velocity

@inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields, calcite_name = Val(:CaCO₃)) = (
  - phytoplankton_primary_production(bgc, i, j, k, fields, auxiliary_fields)
  + plankton_inorganic_carbon_waste(bgc, i, j, k, fields, auxiliary_fields) 
  + detritus_inorganic_carbon_waste(bgc, i, j, k, fields, auxiliary_fields)
  + calcite_dissolution(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
  - calcite_precipitation(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
)

@inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields, calcite_name = Val(:CaCO₃)) = (
    bgc(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields) * (1 - 1/16)
  - bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields) * (1 + 1/16) 
  - 2 * calcite_precipitation(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
  + 2 * calcite_dissolution(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
)

@inline (bgc::NutrientsPlanktonDetritus{<:Nutrient, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields, calcite_name = Val(:CaCO₃)) = (
    bgc(i, j, k, grid, Val(:N), clock, fields, auxiliary_fields) # idk what the correct ratio for the phosphate is
  - 2 * calcite_precipitation(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
  + 2 * calcite_dissolution(bgc, calcite_name, i, j, k, fields, auxiliary_fields)
)

# implicit calcite
@inline calcite_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, calcite_name, args...) = (
    biological_calcite_dissolution(bgc, args...)
  + biological_calcite_loss(bgc, args...) # skips the calcite pool and instantly dissolved
)

@inline calcite_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, calcite_name, args...) = 
    biological_calcite_precipitation(bgc, args...)

# explicit calcite
@inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, val_name::Val{:CaCO₃}, clock, args...) = (
    calcite_production(bgc, val_name, i, j, k, args...)
  - kinetic_dissolution(bgc, val_name, i, j, k, args...)
)

@inline calcite_production(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}}, val_name, args...) = (
    biological_calcite_loss(bgc, args...) 
  + kinetic_precipitation(bgc, val_name, args...)
)

@inline calcite_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}}, val_name, args...) = (
    biological_calcite_dissolution(bgc, args...)
  + kinetic_dissolution(bgc, val_name, args...)
)

@inline calcite_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}}, calcite_name, args...) = (
    biological_calcite_precipitation(bgc, args...)
  + kinetic_precipitation(bgc, calcite_name, args...)
)

@inline function _kinetic_precipitation(bgc, i, j, k, calcite_saturation_state)
    r = bgc.carbonate_system.precipitation_rate_constant
    n = bgc.carbonate_system.precipitation_order

    Ω = @inbounds calcite_saturation_state[i, j, k]

    excess = max(zero(Ω), Ω - one(Ω))
    
    return r * excess ^ n
end

@inline function _kinetic_dissolution(bgc, i, j, k, calcite_saturation_state, calcite)
    r = bgc.carbonate_system.dissolution_rate_constant
    n = bgc.carbonate_system.dissolution_order

    @inbounds begin
        Ω = calcite_saturation_state[i, j, k]
        CaCO₃ = calcite[i, j, k]
    end

    excess = max(zero(Ω), one(Ω) - Ω)
    
    return r * excess ^ n * CaCO₃
end

@inline kinetic_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{1, true}},
                              ::Val{:CaCO₃}, i, j, k, fields, auxiliary_fields) =
    @inbounds _kinetic_precipitation(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[1])

@inline kinetic_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{1, true}},
                            ::Val{:CaCO₃}, i, j, k, fields, auxiliary_fields) =
    @inbounds _kinetic_dissolution(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[1], fields.CaCO₃)

function manifest_carbonate_replicates!(N)
    if N>1
      for n in 1:N
          DIC_name = Symbol(:DIC, n)
          Alk_name = Symbol(:Alk, n)
          CaCO₃_name = Symbol(:CaCO₃, n)
          @eval begin
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:DIC), clock, fields, auxiliary_fields, Val($(QuoteNode(CaCO₃_name))))
              
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:Alk), clock, fields, auxiliary_fields, Val($(QuoteNode(CaCO₃_name))))
              
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, val_name::Val{$(QuoteNode(CaCO₃_name))}, clock, args...) = (
                  calcite_production(bgc, val_name, i, j, k, args...)
                - kinetic_dissolution(bgc, val_name, i, j, k, args...)
              )

              @inline kinetic_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}},
                                            ::Val{$(QuoteNode(CaCO₃_name))},
                                            i, j, k, fields, auxiliary_fields) =
                  @inbounds _kinetic_precipitation(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[$(n)])

              @inline kinetic_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}},
                                          ::Val{$(QuoteNode(CaCO₃_name))},
                                          i, j, k, fields, auxiliary_fields) =
                  @inbounds _kinetic_dissolution(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[$(n)], getproperty(fields, $(QuoteNode(CaCO₃_name))))

              @inline biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, ::Val{$(QuoteNode(CaCO₃_name))}) =
                  biogeochemical_drift_velocity(bgc, Val(:CaCO₃))
          end
      end
    end
end