"""
    Abiotic
No life, does nothing
"""
struct Abiotic end

required_biogeochemical_tracers(::Abiotic) = tuple()
required_biogeochemical_auxiliary_fields(::Abiotic) = tuple()

@inline limiting_nutrients(::Abiotic, args...) = tuple()
@inline inorganic_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline dissolved_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline solid_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)

# thats annoying, need these methods for ambiguity
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:N}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:NO₃}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:NH₄}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:PO₄}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:Fe}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, ::Val{:Si}, args...) where FT = zero(FT)

Base.summary(::Abiotic) = string("Abiotic")
Base.show(io::IO, ab::Abiotic) = print(io, summary(ab))