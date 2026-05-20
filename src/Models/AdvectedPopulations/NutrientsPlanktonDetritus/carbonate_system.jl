using Oceananigans: defaults, Units

using Oceananigans.Architectures: architecture
using Oceananigans.Fields: CenterField, Center
using Oceananigans.Grids: znode, AbstractGrid, Flat
using Oceananigans.Utils: launch!

using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

#using OceanBioME.
using KernelAbstractions: @kernel, @index

import Oceananigans.Biogeochemistry: update_biogeochemical_state!

"""
    CarbonateSystem{N, CAL}

`CarbonateSystem` defines the carbonate system for the `LOBSTER` biogeochemical
model and evolves dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`).
`DIC` is the total inorganic carbon, typically considered to be dissolved 
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

`Alk`alinity is produced (by the removal of nitrose acid (?)), and removed by
the uptake of calcite into phytoplankton.

`DIC` and `Alk` concentration is only one way coupled with the rest of the 
biogeochemistry and *does not* effect any other groups (e.g. acidifcation does
not effect phytoplankton growth). To capture this effect a different `plankton` 
could be defined.

Multiple (N) instances of the carbonate system can evolve in parallel and the tracers
will be named `DIC1` ... `DICN` and `Alk1` ... `AlkN`. This is setup by passing an 
integer argument when the model is constructed like `CarbonateSystem(2)`.

You may get method overwrite warnings if you repeatedly define carbonate systems
with N>1, this shouldn't cause problems.
"""
struct CarbonateSystem{N, CAL, CSS, CC, FT}
      calcite_saturation_state :: CSS
              carbon_chemistry :: CC
    gravitational_acceleration :: FT
     dissolution_rate_constant :: FT 
             dissolution_order :: FT
   precipitation_rate_constant :: FT
           precipitation_order :: FT
end

function CarbonateSystem(replicates = 1; 
                         grid = nothing, 
                         explicit_calcite = false,
                         carbon_chemistry::CC = CarbonChemistry(),
                         gravitational_acceleration = defaults.gravitational_acceleration,
                         dissolution_rate_constant = 1e-4 / day,
                         dissolution_order = 1.0,
                         precipitation_rate_constant = 1e-5 / day,
                         precipitation_order = 2.0) where CC
                         
    manifest_carbonate_replicates!(replicates)

    if explicit_calcite
        calcite_saturation_state = repeat([CenterField(grid)], replicates)
    else
        calcite_saturation_state = nothing
    end

    CSS = typeof(calcite_saturation_state)
    FT = eltype(grid)

    return CarbonateSystem{replicates, explicit_calcite, CSS, CC, FT}(calcite_saturation_state, 
                                                                      carbon_chemistry, 
                                                                      gravitational_acceleration,
                                                                      dissolution_rate_constant,
                                                                      dissolution_order,
                                                                      precipitation_rate_constant,
                                                                      precipitation_order)
end

update_biogeochemical_state!(model, ::CarbonateSystem{<:Any, Nothing}) = nothing

function update_biogeochemical_state!(model, bgc::CarbonateSystem{N}) where N
    carbon_chemistry = bgc.carbon_chemistry
    grid = model.grid

    T = model.tracers.T
    S = model.tracers.S

    arch = architecture(grid)

    for n in 1:N
        DIC = model.tracers[Symbol(:DIC, ifelse(n==1, "", n))]
        Alk = model.tracers[Symbol(:Alk, ifelse(n==1, "", n))]
        Ω   = @inbounds bgc.calcite_saturation_state[n]

        launch!(arch, model.grid, :xyz, 
                compute_calcite_saturation!, 
                Ω, DIC, Alk, T, S, grid, carbon_chemistry, bgc)
    end
end

@inline function pressure(i, j, k, grid, g, ρ)
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

required_biogeochemical_tracers(::CarbonateSystem{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::CarbonateSystem{1, true}) = (:DIC, :Alk, :CaCO₃)

required_biogeochemical_tracers(::CarbonateSystem{N}) where N = (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)
required_biogeochemical_tracers(::CarbonateSystem{N}) where N = (map(n->Symbol(:DIC, n), 1:N)..., 
                                                                 map(n->Symbol(:Alk, n), 1:N)...,
                                                                 map(n->Symbol(:CaCO₃, n), 1:N)...)

@inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields) = (
  - phytoplankton_primary_production(bgc, i, j, k, fields, auxiliary_fields)
  + plankton_inorganic_carbon_waste(bgc, i, j, k, fields, auxiliary_fields)
  + detritus_inorganic_carbon_waste(bgc, i, j, k, fields, auxiliary_fields)
  - biological_calcite_storage(bgc, i, j, k, fields, auxiliary_fields)
  + calcite_dissolution(bgc, i, j, k, fields, auxiliary_fields)
)

@inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields) = (
    bgc(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields) * (1 - 1/16)
  - bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields) * (1 + 1/16) 
  - 2 * biological_calcite_storage(bgc, i, j, k, fields, auxiliary_fields)
  + 2 * calcite_dissolution(bgc, i, j, k, fields, auxiliary_fields)
)

@inline (bgc::NutrientsPlanktonDetritus{<:Nutrient, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields) = (
    bgc(i, j, k, grid, Val(:N), clock, fields, auxiliary_fields)
  - 2 * biological_calcite_storage(bgc, i, j, k, fields, auxiliary_fields)
  + 2 * calcite_dissolution(bgc, i, j, k, fields, auxiliary_fields)
)

@inline calcite_production(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, args...) =
    biological_calcite_production(bgc, args...)

@inline calcite_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem}, args...) =
    biological_calcite_dissolution(bgc, args...)

# explicit calcite
@inline (bgc::NutrientsPlanktonDetritus{<:Nutrient, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:CaCO₃}, args...) = (
    calcite_production(bgc, i, j, k, grid, args...)
  - calcite_dissolution(bgc, i, j, k, grid, args...)
)

@inline calcite_production(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}}, args...) = (
    biological_calcite_production(bgc, args...) 
  + kinetic_precipitation(bgc, args...)
)

@inline calcite_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{<:Any, true}}, args...) = (
    biological_calcite_dissolution(bgc, args...)
  + kinetic_dissolution(bgc, args...)
)

@inline function _kinetic_precipitation(bgc, i, j, k, calcite_saturation_state)
    r = bgc.carbonate_system.precipitation_rate_constant
    n = bgc.carbonate_system.precipitation_order

    Ω = @inbounds calcite_saturation_state[i, j, k]

    excess = max(zero(calcite_saturation_state), Ω - one(calcite_saturation_state))
    
    return r * excess ^ n
end

@inline function _kinetic_dissolution(bgc, i, j, k, calcite_saturation_state, calcite)
    r = bgc.carbonate_system.dissolution_rate_constant
    n = bgc.carbonate_system.dissolution_order

    @inbounds begin
        Ω = calcite_saturation_state[i, j, k]
        CaCO₃ = calcite[i, j, k]
    end

    excess = max(zero(calcite_saturation_state), one(calcite_saturation_state) - Ω)
    
    return r * excess ^ n * CaCO₃
end

@inline kinetic_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{1, true}},
                              i, j, k, fields, auxiliary_fields) =
    @inbounds _kinetic_precipitation(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[1])

@inline kinetic_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{1, true}},
                            i, j, k, fields, auxiliary_fields) =
    @inbounds _kinetic_dissolution(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[1], fields.CaCO₃)

function manifest_carbonate_replicates!(N)
    if N>1
      for n in 1:N
          DIC_name = Symbol(:DIC, n)
          Alk_name = Symbol(:Alk, n)
          CaCO₃_name = Symbol(:CaCO₃, n)
          @eval begin
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:DIC), clock, fields, auxiliary_fields)
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:Alk), clock, fields, auxiliary_fields)
              @inline (bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{$(QuoteNode(CaCO₃_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:CaCO₃), clock, fields, auxiliary_fields)

              @inline kinetic_precipitation(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{$(QuoteNode(n)), true}},
                                            i, j, k, fields, auxiliary_fields) =
                  @inbounds _kinetic_precipitation(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[$(QuoteNode(n))])

              @inline kinetic_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonateSystem{$(QuoteNode(n)), true}},
                                          i, j, k, fields, auxiliary_fields) =
                  @inbounds _kinetic_dissolution(bgc, i, j, k, bgc.carbonate_system.calcite_saturation_state[$(QuoteNode(n))], fields[$(QuoteNode(CaCO₃_name))])
          end
      end
    end
end