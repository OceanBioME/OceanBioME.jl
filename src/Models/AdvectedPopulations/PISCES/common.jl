using KernelAbstractions: @kernel, @index

using Oceananigans.Grids: znode, zspacing

@inline shear(z, zₘₓₗ, background_shear, mixed_layer_shear) = ifelse(z <= zₘₓₗ, background_shear, mixed_layer_shear) # Given as 1 in Aumont paper

@inline latitude(φ, y) = φ
@inline latitude(::Nothing, y) = y

# we should probably extend this to use DateTime dates at some point
@inline function day_length(φ, t)
    # as per Forsythe et al., 1995 (https://doi.org/10.1016/0304-3800(94)00034-F)
    p = asind(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (floor(Int, t / day) - 186)))))

    return 24 - 24 / 180 * acosd((sind(0.8333) + sind(φ) * sind(p)) / (cosd(φ) * cosd(p)))
end

@kwdef struct DepthDependantSinkingSpeed{FT}
    minimum_speed :: FT = 2/day
    maximum_speed :: FT = 200/day
    maximum_depth :: FT = 5000.0
end

# I can't find any explanation as to why this might depend on the euphotic depth
function (p::DepthDependantSinkingSpeed)(i, j, k, grid, mixed_layer_depth, euphotic_depth)
    zₘₓₗ = @inbounds mixed_layer_depth[i, j, k]
    zₑᵤ  = @inbounds euphotic_depth[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    return - p.minimum_speed + (p.maximum_speed - p.minimum_speed) * min(0, z - min(zₘₓₗ, zₑᵤ)) / 5000
end

compute_mean_mixed_layer_vertical_diffusivity!(κ::ConstantField, mixed_layer_depth, model) = nothing

compute_mean_mixed_layer_vertical_diffusivity!(κ, mixed_layer_depth, model) =
    compute_mean_mixed_layer_vertical_diffusivity!(model.closure, κ, mixed_layer_depth, model.diffusivity_fields)

# if no closure is defined we just assume its pre-set 
compute_mean_mixed_layer_vertical_diffusivity!(closure::Nothing, 
                                               mean_mixed_layer_vertical_diffusivity, 
                                               mixed_layer_depth, 
                                               diffusivity_fields) = nothing

function compute_mean_mixed_layer_vertical_diffusivity!(closure, mean_mixed_layer_vertical_diffusivity, mixed_layer_depth, diffusivity_fields)
    # this is going to get messy
    κ = phytoplankton_diffusivity(closure, diffusivity_fields)

    launch!(arch, grid, :xy, _compute_mean_mixed_layer_vertical_diffusivity!, mean_mixed_layer_vertical_diffusivity, mixed_layer_depth, κ)

    fill_halo_regions!(mean_mixed_layer_vertical_diffusivity)

    return nothing
end

@kernel function _compute_mean_mixed_layer_vertical_diffusivity!(κₘₓₗ, mixed_layer_depth, κ)
    i, j = @index(Global, NTuple)

    zₘₓₗ = @inbounds mixed_layer_depth[i, j]

    @inbounds κₘₓₗ[i, j] = 0

    integration_depth = 0

    for k in grid.Nz:-1:1
        if znode(i, j, k, grid, Center(), Center(), Center()) > zₘₓₗ
            Δz = zspacing(i, j, k, grid, Center(), Center(), Center())
            κₘₓₗ[i, j] += κ[i, j, k] * Δz # I think sometimes vertical diffusivity is face located?
            integration_depth += Δz
        end
    end

    κₘₓₗ[i, j] /= integration_depth
end

# this does not belong here - lets add them when a particular closure is needed
using Oceananigans.TurbulenceClosures: ScalarDiffusivity, ScalarBiharmonicDiffusivity, VerticalFormulation, ThreeDimensionalFormulation, formulation

phytoplankton_diffusivity(closure, diffusivity_fields) =
    phytoplankton_diffusivity(formulation(closure), closure, diffusivity_fields)

phytoplankton_diffusivity(closure::Tuple, diffusivity_fields) = 
    sum(map(n -> phytoplankton_diffusivity(closure[n], diffusivity_fields[n]), 1:length(closure)))

phytoplankton_diffusivity(formulation, closure, diffusivit_fields) = ZeroField()

const NotHorizontalFormulation = Union{VerticalFormulation, ThreeDimensionalFormulation}

phytoplankton_diffusivity(::NotHorizontalFormulation, closure, diffusivity_fields) = 
    throw(ErrorException("Mean mixed layer vertical diffusivity can not be calculated for $(closure)"))

phytoplankton_diffusivity(::NotHorizontalFormulation, 
                          closure::Union{ScalarDiffusivity, ScalarBiharmonicDiffusivity}, 
                          diffusivity_fields) = 
    phytoplankton_diffusivity(closure.κ)

phytoplankton_diffusivity(diffusivity_field) = diffusivity_field
phytoplankton_diffusivity(diffusivity_field::Number) = ConstantField(diffusivity_field)
phytoplankton_diffusivity(diffusivity_fields::NamedTuple) = phytoplankton_diffusivity(diffusivity_fields.P)
phytoplankton_diffusivity(::Function) = 
    throw(ErrorException("Can not compute mean mixed layer vertical diffusivity for `Function` type diffusivity, changing to a `FunctionField` would work"))


