#####
##### Mean mixed layer diffusivity
#####
compute_mean_mixed_layer_vertical_diffusivity!(κ::ConstantField, mixed_layer_depth, model) = nothing

compute_mean_mixed_layer_vertical_diffusivity!(κ, mixed_layer_depth, model) =
    compute_mean_mixed_layer_vertical_diffusivity!(model.closure, κ, mixed_layer_depth, model.diffusivity_fields, grid)

# if no closure is defined we just assume its pre-set 
compute_mean_mixed_layer_vertical_diffusivity!(closure::Nothing, 
                                               mean_mixed_layer_vertical_diffusivity, 
                                               mixed_layer_depth, 
                                               diffusivity_fields, grid) = nothing

function compute_mean_mixed_layer_vertical_diffusivity!(closure, mean_mixed_layer_vertical_diffusivity, mixed_layer_depth, diffusivity_fields, grid)
    # this is going to get messy
    κ = phytoplankton_diffusivity(closure, diffusivity_fields)

    launch!(arch, grid, :xy, _compute_mean_mixed_layer_vertical_diffusivity!, mean_mixed_layer_vertical_diffusivity, mixed_layer_depth, κ, grid)

    fill_halo_regions!(mean_mixed_layer_vertical_diffusivity)

    return nothing
end

@kernel function _compute_mean_mixed_layer_vertical_diffusivity!(κₘₓₗ, mixed_layer_depth, κ, grid)
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

#####
##### Mean mixed layer light
#####

function compute_mean_mixed_layer_light!(mean_PAR, mixed_layer_depth, PAR, model)
    grid = model.grid
    
    launch!(arch, grid, :xy, _compute_mean_mixed_layer_light!, mean_PAR, mixed_layer_depth, PAR, grid)

    fill_halo_regions!(mean_PAR)

    return nothing
end

@kernel function _compute_mean_mixed_layer_light!(mean_PAR, mixed_layer_depth, PAR, grid)
    i, j = @index(Global, NTuple)

    zₘₓₗ = @inbounds mixed_layer_depth[i, j]

    @inbounds mean_PAR[i, j] = 0

    integration_depth = 0

    for k in grid.Nz:-1:1
        if znode(i, j, k, grid, Center(), Center(), Center()) > zₘₓₗ
            Δz = zspacing(i, j, k, grid, Center(), Center(), Center())
            mean_PAR[i, j] += PAR[i, j, k] * Δz 
            integration_depth += Δz
        end
    end

    mean_PAR[i, j] /= integration_depth
end

#####
##### Informaiton about diffusivity fields
#####

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


