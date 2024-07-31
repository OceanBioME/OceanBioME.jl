using Oceananigans.Utils: SumOfArrays

include("morel_coefficients.jl")

struct ScaledSurfaceFunction{SF, SC} <: Function
    surface_function :: SF
        scale_factor :: SC
end

@inline (ssf::ScaledSurfaceFunction)(X...) = ssf.surface_function(X...) * ssf.scale_factor

struct MultiBandPhotosyntheticallyActiveRadiation{F, FN, K, E, C, SPAR}
    fields :: F
    field_names :: FN

    water_attenuation_coefficient :: K
    chlorophyll_exponent :: E
    chlorophyll_attenuation_coefficient :: C

    surface_PAR :: SPAR
end

"""
    MultiBandPhotosyntheticallyActiveRadiation(; )

Keyword Arguments
==================

- `grid`: grid for building the model on
- `water_red_attenuation`, ..., `phytoplankton_chlorophyll_ratio`: parameter values
- `surface_PAR`: function (or array in the future) for the photosynthetically available radiation at the surface, 
   which should be `f(x, y, t)` where `x` and `y` are the native coordinates (i.e. meters for rectilinear grids
   and latitude/longitude as appropriate)
"""
function MultiBandPhotosyntheticallyActiveRadiation(; grid, 
                                                      bands = ((400, 500), (500, 600), (600, 700)), #nm
                                                      base_bands = MOREL_λ,
                                                      base_water_attenuation_coefficient = MOREL_kʷ,
                                                      base_chlorophyll_exponent = MOREL_e,
                                                      base_chlorophyll_attenuation_coefficient = MOREL_χ,
                                                      field_names = [par_symbol(n, length(bands)) for n in 1:length(bands)],
                                                      surface_PAR = default_surface_PAR)
    Nbands = length(bands)

    kʷ = zeros(eltype(grid), Nbands)
    e  = zeros(eltype(grid), Nbands)
    χ  = zeros(eltype(grid), Nbands)

    for (n, band) in enumerate(bands)
        idx1 = findlast(base_bands .<= band[1])
        idx2 = findlast(base_bands .< band[2])

        kʷ[n] = sum(base_water_attenuation_coefficient[idx1:idx2] .* base_bands[idx1:idx2]) / sum(base_bands[idx1:idx2])
        e[n]  = sum(base_chlorophyll_exponent[idx1:idx2] .* base_bands[idx1:idx2]) / sum(base_bands[idx1:idx2])
        χ[n]  = sum(base_chlorophyll_attenuation_coefficient[idx1:idx2] .* base_bands[idx1:idx2]) / sum(base_bands[idx1:idx2])
    end
    

    wrapped_surface_PAR_function = ScaledSurfaceFunction(surface_PAR, 1 / Nbands)

    fields = [CenterField(grid; boundary_conditions = 
                            regularize_field_boundary_conditions(
                                FieldBoundaryConditions(top = ValueBoundaryCondition(wrapped_surface_PAR_function)), grid, name)) for name in field_names]


    return MultiBandPhotosyntheticallyActiveRadiation(fields, field_names, kʷ, e, χ, surface_PAR)
end

function par_symbol(n, nbands)
    if nbands <= 9
        return Symbol(:PAR, Char('\xe2\x82\x80'+n))
    else
        Symbol(:PAR, n)
    end
end

@kernel function update_MultiBandPhotosyntheticallyActiveRadiation!(PAR_model, grid, Chl, t) 
    i, j = @index(Global, NTuple)

    k, k′ = domain_boundary_indices(RightBoundary(), grid.Nz)

    X = z_boundary_node(i, j, k′, grid, Center(), Center())

    surface_PAR = PAR_model.surface_PAR(X..., t)
    fields = PAR_model.fields

    kʷ = PAR_model.water_attenuation_coefficient
    e  = PAR_model.chlorophyll_exponent
    χ  = PAR_model.chlorophyll_attenuation_coefficient

    Nbands = length(fields)

    zᶜ = znodes(grid, Center(), Center(), Center())

    # first point below surface
    k = grid.Nz
    for n in 1:Nbands
        @inbounds fields[n][i, j, k] = @inbounds (surface_PAR / Nbands) * exp(zᶜ[grid.Nz] * (kʷ[n] + χ[n] * Chl[i, j, k] ^ e[n]))
    end

    # the rest of the points
    for k in grid.Nz-1:-1:1
        Δz = @inbounds zᶜ[k] - zᶜ[k + 1] 
        for n in 1:Nbands
            @inbounds fields[n][i, j, k] = @inbounds fields[n][i, j, k + 1] * exp(Δz * (kʷ[n] + χ[n] * Chl[i, j, k] ^ e[n]))
        end
    end
end


function update_biogeochemical_state!(model, PAR::MultiBandPhotosyntheticallyActiveRadiation)
    grid = model.grid

    arch = architecture(grid)

    launch!(arch, grid, :xy, update_MultiBandPhotosyntheticallyActiveRadiation!, PAR, grid, chlorophyll(model.biogeochemistry, model), model.clock.time)

    for field in PAR.fields
        fill_halo_regions!(field, model.clock, fields(model))
    end
end

summary(par::MultiBandPhotosyntheticallyActiveRadiation) = 
    string("Multi band light attenuation model with $(length(par.fields)) bands $(par.field_names)")
show(io::IO, model::MultiBandPhotosyntheticallyActiveRadiation) = print(io, summary(model))

biogeochemical_auxiliary_fields(par::MultiBandPhotosyntheticallyActiveRadiation) = 
    NamedTuple{(:PAR, par.field_names...)}([SumOfArrays{length(par.fields)}(par.fields...), par.fields...])

adapt_structure(to, par::MultiBandPhotosyntheticallyActiveRadiation) = 
    MultiBandPhotosyntheticallyActiveRadiation(adapt(to, par.fields),
                                               nothing, # don't need this in the kernel
                                               adapt(to, par.water_attenuation_coefficient),
                                               adapt(to, chlorophyll_exponent),
                                               adapt(to, chlorophyll_attenuation_coefficient),
                                               adapt(to, surface_PAR))
