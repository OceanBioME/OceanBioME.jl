using Oceananigans.Fields: Field, OneField, CenterField
using Oceananigans.Grids: Center
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Architectures: architecture
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields
using Oceananigans.Utils: launch!
using Oceananigans: fields
using KernelAbstractions: @kernel, @index

using Adapt

floor_index_field(grid) = OneField(Int)

function floor_index_field(grid::ImmersedBoundaryGrid)
    floor_indices = Field{Center, Center, Nothing}(grid, Int)
    launch!(architecture(grid), grid, :xy, _compute_floor_indices!, grid, floor_indices, size(grid, 3))
    return floor_indices
end

@kernel function _compute_floor_indices!(grid, floor_indices, Nz)
    i, j = @index(Global, NTuple)

    k = 1
    while immersed_cell(i, j, k, grid.underlying_grid, grid.immersed_boundary) && k < Nz
        k += 1
    end

    @inbounds floor_indices[i, j, 1] = k
end

#####
##### Sinking types
#####

struct ExplicitSinking{SS}
    sinking_speeds :: SS
end

Adapt.adapt_structure(to, s::ExplicitSinking) =
    ExplicitSinking(Adapt.adapt(to, s.sinking_speeds))

struct ImplicitSinking{FT, RM, FL, FI}
    dissolution_length :: FT
    remineralisation    :: RM
    floor_flux          :: FL
    floor_indices       :: FI
    open_bottom         :: Bool
end

function ImplicitSinking(grid, dissolution_length::Number; tracer_names, open_bottom = true)
    FT = eltype(grid)
    fi = floor_index_field(grid)
    remin = NamedTuple{tracer_names}(ntuple(length(tracer_names)) do _
        f = CenterField(grid)
        fill!(f, 0)
        f
    end)
    ff = NamedTuple{tracer_names}(ntuple(length(tracer_names)) do _
        f = Field{Center, Center, Nothing}(grid)
        fill!(f, 0)
        f
    end)
    return ImplicitSinking(convert(FT, dissolution_length), remin, ff, fi, open_bottom)
end

function ImplicitSinking(grid, dissolution_length::NamedTuple; open_bottom = true)
    FT = eltype(grid)
    tracer_names = keys(dissolution_length)
    fi = floor_index_field(grid)
    dl = NamedTuple{tracer_names}(convert.(FT, values(dissolution_length)))
    remin = NamedTuple{tracer_names}(ntuple(length(tracer_names)) do _
        f = CenterField(grid)
        fill!(f, 0)
        f
    end)
    ff = NamedTuple{tracer_names}(ntuple(length(tracer_names)) do _
        f = Field{Center, Center, Nothing}(grid)
        fill!(f, 0)
        f
    end)
    return ImplicitSinking(dl, remin, ff, fi, open_bottom)
end

Adapt.adapt_structure(to, s::ImplicitSinking) =
    ImplicitSinking(s.dissolution_length,
                    Adapt.adapt(to, s.remineralisation),
                    Adapt.adapt(to, s.floor_flux),
                    Adapt.adapt(to, s.floor_indices),
                    s.open_bottom)

Base.summary(::ExplicitSinking) = "ExplicitSinking"
Base.summary(s::ImplicitSinking) = "ImplicitSinking(ℓ=$(s.dissolution_length) m)"

#####
##### Shared implicit-sinking column sweep
#####

"""
    implicit_sinking_production(i, j, k, grid, detritus, bgc, fields, aux, val_name)

Return the source of sinking particulate matter at cell `(i,j,k)` for the tracer
identified by `val_name` (a `Val{:name}`). Each detritus type adds methods for its
particulate tracers.
"""
function implicit_sinking_production end

dissolution_length(s::ImplicitSinking{<:Number}, name) = s.dissolution_length
dissolution_length(s::ImplicitSinking, name) = s.dissolution_length[name]

@kernel function implicit_sinking_column!(grid, detritus, bgc, model_fields, aux,
                                           remin, floor_flux, floor_indices,
                                           ℓ, open_bottom, Nz, val_name)
    i, j = @index(Global, NTuple)

    FT = eltype(grid)
    kf = @inbounds floor_indices[i, j, 1]
    F = zero(FT)

    for k in Nz:-1:kf
        Δz = Δzᶜᶜᶜ(i, j, k, grid)

        Π = implicit_sinking_production(i, j, k, grid, detritus, bgc, model_fields, aux, val_name)

        F_in = F
        F = (F_in + Π * Δz) / (one(FT) + Δz / ℓ)

        at_closed_floor = (k == kf) & !open_bottom
        R = ifelse(at_closed_floor, Π + F_in / Δz, Π + (F_in - F) / Δz)
        F = ifelse(at_closed_floor, zero(FT), F)

        @inbounds remin[i, j, k] = R
    end

    @inbounds floor_flux[i, j, 1] = ifelse(open_bottom, F, zero(FT))
end

