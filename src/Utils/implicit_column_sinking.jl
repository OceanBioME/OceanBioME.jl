using Oceananigans.Fields: Field, OneField, CenterField
using Oceananigans.Grids: Center, znode
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

using Adapt

floor_index_field(grid) = OneField(Int)

function floor_index_field(grid::ImmersedBoundaryGrid)
    floor_indices = Field{Center, Center, Nothing}(grid, Int)
    launch!(architecture(grid), grid, :xy, compute_floor_indices!, grid, floor_indices, size(grid, 3))
    return floor_indices
end

@kernel function compute_floor_indices!(grid, floor_indices, Nz)
    i, j = @index(Global, NTuple)

    k = 1
    while immersed_cell(i, j, k, grid.underlying_grid, grid.immersed_boundary) && k < Nz
        k += 1
    end

    @inbounds floor_indices[i, j, 1] = k
end

"""
    ExplicitSinking(sinking_speeds)

Sinking handled by the tracer advection scheme via `biogeochemical_drift_velocity`.
`sinking_speeds` is a `NamedTuple` of velocity fields `(u, v, w)`.
"""
struct ExplicitSinking{SS}
    sinking_speeds :: SS
end

Adapt.adapt_structure(to, s::ExplicitSinking) =
    ExplicitSinking(Adapt.adapt(to, s.sinking_speeds))

"""
    ImplicitSinking(grid, dissolution_length; tracer_names = keys(dissolution_length), open_bottom = true)

Sinking handled as an implicit vertical redistribution with an exponential profile characterised
by `dissolution_length` (in metres). This avoids the CFL restriction of explicit advection for
fast-sinking particles. When `open_bottom = true` material reaching the bottom cell leaves the
domain; otherwise it accumulates.
"""
struct ImplicitSinking{FT, RM, FL, FI}
    dissolution_length :: FT
      remineralisation :: RM
            floor_flux :: FL
         floor_indices :: FI
           open_bottom :: Bool
end

convert_dissolution_length(FT, ℓ::Number) = convert(FT, ℓ)
convert_dissolution_length(FT, ℓ::NamedTuple) = map(ℓ -> convert_dissolution_length(FT, ℓ), ℓ)
convert_dissolution_length(FT, ℓ) = ℓ

function ImplicitSinking(grid, dissolution_length; tracer_names = keys(dissolution_length), open_bottom = true)
    remineralisation = NamedTuple{tracer_names}(map(_ -> CenterField(grid), tracer_names))
    floor_flux = NamedTuple{tracer_names}(map(_ -> Field{Center, Center, Nothing}(grid), tracer_names))

    return ImplicitSinking(convert_dissolution_length(eltype(grid), dissolution_length),
                           remineralisation,
                           floor_flux,
                           floor_index_field(grid),
                           open_bottom)
end

Adapt.adapt_structure(to, s::ImplicitSinking) =
    ImplicitSinking(Adapt.adapt(to, s.dissolution_length),
                    Adapt.adapt(to, s.remineralisation),
                    Adapt.adapt(to, s.floor_flux),
                    Adapt.adapt(to, s.floor_indices),
                    s.open_bottom)

Base.summary(::ExplicitSinking) = "ExplicitSinking"
Base.summary(s::ImplicitSinking) = "ImplicitSinking(ℓ=$(s.dissolution_length) m)"

"""
    DepthDependentDissolutionLength(f)

Wrap a function of depth `z -> ℓ(z)` so that it can be used as a spatially varying dissolution
length in [`ImplicitSinking`](@ref). The wrapped function receives the cell-centre depth (negative
below the surface) and must return the dissolution length in metres.

Example
=======

```julia
Detritus(grid; dissolution_length = DepthDependentDissolutionLength(z -> 100 + 2 * abs(z)))
```
"""
struct DepthDependentDissolutionLength{F}
    f :: F
end

@inline (d::DepthDependentDissolutionLength)(i, j, k, grid, clock, fields) =
    d.f(znode(i, j, k, grid, Center(), Center(), Center()))

Adapt.adapt_structure(to, d::DepthDependentDissolutionLength) =
    DepthDependentDissolutionLength(Adapt.adapt(to, d.f))

@inline dissolution_length(ℓ::Number, i, j, k, grid, clock, fields) = ℓ
@inline dissolution_length(ℓ, i, j, k, grid, clock, fields) = ℓ(i, j, k, grid, clock, fields)

"""
    implicit_sinking_production(i, j, k, grid, detritus, bgc, fields, aux, val_name)

Return the source of sinking particulate matter at cell `(i,j,k)` for the tracer
identified by `val_name` (a `Val{:name}`). Each detritus type adds methods for its
particulate tracers.
"""
function implicit_sinking_production end

dissolution_length(s::ImplicitSinking{<:NamedTuple}, name) = s.dissolution_length[name]
dissolution_length(s::ImplicitSinking, name) = s.dissolution_length

@kernel function implicit_sinking_column!(grid, detritus, bgc, model_fields, aux,
                                           remineralisation, floor_flux, floor_indices,
                                           ℓ, open_bottom, Nz, val_name, clock)
    i, j = @index(Global, NTuple)

    FT = eltype(grid)
    kf = @inbounds floor_indices[i, j, 1]
    F = zero(FT)

    for k in Nz:-1:kf
        Δz = Δzᶜᶜᶜ(i, j, k, grid)

        Π = implicit_sinking_production(i, j, k, grid, detritus, bgc, model_fields, aux, val_name)

        ℓₖ = dissolution_length(ℓ, i, j, k, grid, clock, model_fields)

        F_in = F
        F = (F_in + Π * Δz) / (one(FT) + Δz / ℓₖ)

        at_closed_floor = (k == kf) & !open_bottom
        R = ifelse(at_closed_floor, Π + F_in / Δz, Π + (F_in - F) / Δz)
        F = ifelse(at_closed_floor, zero(FT), F)

        @inbounds remineralisation[i, j, k] = R
    end

    @inbounds floor_flux[i, j, 1] = ifelse(open_bottom, F, zero(FT))
end
