using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ZFaceField, AbstractField, location, Center, Face, compute!, ZeroField
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Grids: AbstractGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Utils: launch!

import Adapt: adapt_structure, adapt

# nothings for constant fields
const valid_sinking_velocity_locations = ((Center, Center, Face), (Nothing, Nothing, Face), (Nothing, Nothing, Nothing))

function setup_velocity_fields(drift_speeds, grid::AbstractGrid, open_bottom; smoothing_distance = 2, three_D = false)
    drift_velocities = []
    for (name, w) in pairs(drift_speeds)
        if isa(w, Number)
            w_field = ZFaceField(grid)
            profile = [- w * ifelse(open_bottom, 1.0, (1 - exp((1-k) / smoothing_distance))) for k=1:grid.Nz]
            for k=1:grid.Nz
                @inbounds w_field[:, :, k] .= profile[k]
            end
            open_bottom || move_taper_to_seafloor!(w_field, grid, profile)
        elseif isa(w, AbstractField)
            location(w) in valid_sinking_velocity_locations ||
                @warn "The location of the sinking velocity field provided for $name is incorrect, it should be (Center, Center, Face)"

            open_bottom || @warn "The sinking velocity provided for $name is a field and therefore `open_bottom=false` can't be enforced automatically"

            compute!(w)
            w_field = w
        else
            @warn "Sinking speed provided for $name was not a number or field so may be unsiutable"
            w_field = w
        end

        push!(drift_velocities, three_D ? (u = ZeroField(), v = ZeroField(), w = w_field) : w_field)
    end

    return NamedTuple{keys(drift_speeds)}(drift_velocities)
end

move_taper_to_seafloor!(w_field, grid, profile) = nothing

function move_taper_to_seafloor!(w_field, grid::ImmersedBoundaryGrid, profile)
    arch = architecture(grid)
    profile = on_architecture(arch, convert.(eltype(w_field), profile))

    launch!(arch, grid, :xy, _move_taper_to_seafloor!, w_field, grid, profile)

    fill_halo_regions!(w_field)

    return nothing
end

@kernel function _move_taper_to_seafloor!(w, grid, profile)
    i, j = @index(Global, NTuple)

    kᵇ = 1
    while (kᵇ ≤ grid.Nz) && immersed_cell(i, j, kᵇ, grid)
        kᵇ += 1
    end

    for k in 1:grid.Nz
        wₖ = @inbounds profile[max(1, k - kᵇ + 1)]
        @inbounds w[i, j, k] = ifelse(k < kᵇ, zero(wₖ), wₖ)
    end
end

adapt_structure(to, velocities::NamedTuple{(:u, :v, :w), Tuple{AbstractField, AbstractField, AbstractField}}) = NamedTuple{(:u, :v, :w)}(adapt.(to, values(velocities)))

function show_sinking_velocities(sinking_velocities::NamedTuple{T, V}) where {T, V}
    str = ""
    if length(T) == 0
        return str
    elseif length(T) == 1
        return "    └── $(T[1]): $(maximum_sinking(sinking_velocities[1])) to $(minimum_sinking(sinking_velocities[1])) m/s"
    else
        for idx in 1:length(T) - 1
            str *= "    ├── $(T[idx]): $(maximum_sinking(sinking_velocities[idx])) to $(minimum_sinking(sinking_velocities[idx])) m/s \n"
        end
        str *= "    └── $(T[end]): $(maximum_sinking(sinking_velocities[end])) to $(minimum_sinking(sinking_velocities[end])) m/s"

        return str
    end
end

maximum_sinking(velocity) = maximum(velocity)
minimum_sinking(velocity) = minimum(velocity)

maximum_sinking(velocity::NamedTuple) = maximum(velocity.w)
minimum_sinking(velocity::NamedTuple) = minimum(velocity.w)
