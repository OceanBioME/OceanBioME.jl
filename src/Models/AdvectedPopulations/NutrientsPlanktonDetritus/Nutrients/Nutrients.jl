module NutrientsModels 

export Nutrients, SingleTracerNutrient, NitrateAmmonia, N, PO₄, Fe, Si

using Adapt

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio

using ..NutrientsPlanktonDetritusModels:
    inorganic_waste,
    nutrient_uptake,
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import Adapt: adapt_structure
import Base: show, summary

"""
    Nutrients(; nitrogen = nothing, phosphate = nothing, iron = nothing, silicate = nothing)

Group the inorganic nutrient pools that limit plankton growth in a
[`NutrientsPlanktonDetritus`](@ref) model. Each of the four slots may be `nothing` (the nutrient is
not tracked and is implicitly conserved), a [`SingleTracerNutrient`](@ref) (`N`, `PO₄`, `Fe`, or
`Si`), or a more detailed component such as [`NitrateAmmonia`](@ref) for the nitrogen slot.

The tracers added to the model are the union of those required by the non-`nothing` slots.

Keyword Arguments
=================

- `nitrogen`: the nitrogen pool, e.g. `N` (single tracer) or [`NitrateAmmonia`](@ref) (`NO₃`/`NH₄`)
- `phosphate`: the phosphate pool, e.g. `PO₄`
- `iron`: the iron pool, e.g. `Fe`
- `silicate`: the silicate pool, e.g. `Si`
"""
@kwdef struct Nutrients{N, P, F, S}
      nitrogen :: N = nothing
     phosphate :: P = nothing
          iron :: F = nothing
      silicate :: S = nothing
end

"""
    SingleTracerNutrient

An `@enum` of the single-tracer inorganic nutrients that can fill a [`Nutrients`](@ref) slot: `N`
(nitrogen), `PO₄` (phosphate), `Fe` (iron), and `Si` (silicate). Each introduces one tracer of the
same name and is the simplest way to represent a limiting nutrient.
"""
@enum SingleTracerNutrient N PO₄ Fe Si

required_biogeochemical_tracers(nutrients::Nutrients) =
    (required_biogeochemical_tracers(nutrients.nitrogen)...,
     required_biogeochemical_tracers(nutrients.phosphate)...,
     required_biogeochemical_tracers(nutrients.iron)...,
     required_biogeochemical_tracers(nutrients.silicate)...)

required_biogeochemical_tracers(tracer::SingleTracerNutrient) =
    (Symbol(tracer), )

required_biogeochemical_auxiliary_fields(::Nutrients) = tuple()

const NutrientsNPD{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients}

include("prognosis.jl")
include("nitrate_ammonia.jl")
include("defaults.jl")
include("adapt_show_methods.jl")

end # module