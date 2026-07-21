"""
    Abiotic()

An empty plankton component for the `plankton` slot of a [`NutrientsPlanktonDetritus`](@ref) model.
It adds no tracers and produces no biological source or sink terms — all uptake and waste fluxes are
zero. It is the default plankton and is useful for purely abiotic experiments (e.g. carbonate
chemistry or gas exchange with no biology).
"""
struct Abiotic <: AbstractPlankton{tuple()} end

required_biogeochemical_tracers(::Abiotic) = tuple()
required_biogeochemical_auxiliary_fields(::Abiotic) = tuple()

@inline limiting_nutrients(::Abiotic, args...) = tuple()
@inline inorganic_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline dissolved_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline solid_waste(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)

# annoyingly need these methods for ambiguity
@inline nutrient_uptake(i, j, k, grid, ::Val{:N}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Val{:NO₃}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Val{:NH₄}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Val{:PO₄}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Val{:Fe}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)
@inline nutrient_uptake(i, j, k, grid, ::Val{:Si}, ::Abiotic, ::NPD{FT}, args...) where FT = zero(FT)

Base.summary(::Abiotic) = string("Abiotic")
Base.show(io::IO, ab::Abiotic) = print(io, summary(ab))