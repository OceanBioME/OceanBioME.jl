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

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import Adapt: adapt_structure
import Base: show, summary

@kwdef struct Nutrients{N, P, F, S} 
      nitrogen :: N = nothing
     phosphate :: P = nothing
          iron :: F = nothing
      silicate :: S = nothing
end

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