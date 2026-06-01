"""
    Abiotic
No life, does nothing
"""
struct Abiotic end

required_biogeochemical_tracers(::Abiotic) = tuple()
required_biogeochemical_auxiliary_fields(::Abiotic) = tuple()

@inline limiting_nutrients(::Abiotic, args...) = tuple()
@inline inorganic_waste(::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline dissolved_waste(::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline solid_waste(::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)

# thats annoying, need these methods for ambiguity
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:N}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:NO₃}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:NH₄}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:PO₄}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:Fe}, args...) where FT = zero(FT)
@inline nutrient_uptake(::Abiotic, ::NPD{FT}, i, j, k, ::Val{:Si}, args...) where FT = zero(FT)

Base.summary(::Abiotic) = string("Abiotic (no life)")
Base.show(io::IO, ab::Abiotic) = print(io, summary(ab))