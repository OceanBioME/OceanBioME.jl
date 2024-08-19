using Oceananigans.Fields: Field
using Oceananigans.Utils: SumOfArrays

include("morel_coefficients.jl")

struct ScaledSurfaceFunction{SF, SC} <: Function
    surface_function :: SF
        scale_factor :: SC
end

@inline (ssf::ScaledSurfaceFunction)(X...) = ssf.surface_function(X...) * ssf.scale_factor

"""
    MultiBandPhotosyntheticallyActiveRadiation{F, FN, K, E, C, SPAR, SPARD}

Light attenuation model with multiple wave length bands where each band (i) is attenuated like:

    ∂PARᵢ(z)/∂z = PARᵢ(kʷ(i) + χ(i)Chl(z)ᵉ⁽ⁱ⁾) 

Where kʷ(i) is the band specific water attenuation coefficient, e(i) the chlorophyll exponent,
and χ(i) the chlorophyll attenuation coefficient.

When the fields are called with `biogeochemical_auxiliary_fields` an additional field named `PAR`
is also returned which is a sum of the bands.
"""
struct MultiBandPhotosyntheticallyActiveRadiation{T, F, FN, K, E, C, SPAR, SPARD}
                                total :: T
                               fields :: F
                          field_names :: FN

        water_attenuation_coefficient :: K
                 chlorophyll_exponent :: E
  chlorophyll_attenuation_coefficient :: C

                          surface_PAR :: SPAR
                 surface_PAR_division :: SPARD
end

"""
    MultiBandPhotosyntheticallyActiveRadiation(; grid, 
                                                 bands = ((400, 500), (500, 600), (600, 700)), #nm
                                                 base_bands = MOREL_λ,
                                                 base_water_attenuation_coefficient = MOREL_kʷ,
                                                 base_chlorophyll_exponent = MOREL_e,
                                                 base_chlorophyll_attenuation_coefficient = MOREL_χ,
                                                 field_names = [par_symbol(n, length(bands)) for n in 1:length(bands)],
                                                 surface_PAR = default_surface_PAR)

Returns a `MultiBandPhotosyntheticallyActiveRadiation` attenuation model of `surface_PAR` in divided 
into `bands` by `surface_PAR_division`. 

The attenuation morel_coefficients are computed from `base_water_attenuation_coefficient`,
`base_chlorophyll_exponent`, and `base_chlorophyll_attenuation_coefficient` which should be 
arrays of the coefficients at `base_bands` wavelengths. 

The returned `field_names` default to `PAR₁`, `PAR₂`, etc., but may be specified by the user instead.

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
                                                      field_names = [par_symbol(n) for n in 1:length(bands)],
                                                      surface_PAR = default_surface_PAR,
                                                      surface_PAR_division = fill(1 / length(bands), length(bands)))
    Nbands = length(bands)

    kʷ = zeros(eltype(grid), Nbands)
    e  = zeros(eltype(grid), Nbands)
    χ  = zeros(eltype(grid), Nbands)

    for (n, band) in enumerate(bands)
        idx1 = findlast(base_bands .<= band[1])
        idx2 = findlast(base_bands .<= band[2])

        kʷ[n] = numerical_mean(base_bands, base_water_attenuation_coefficient, idx1, idx2)
         e[n] = numerical_mean(base_bands, base_chlorophyll_exponent, idx1, idx2)
         χ[n] = numerical_mean(base_bands, base_chlorophyll_attenuation_coefficient, idx1, idx2)
    end

    sum(surface_PAR_division) == 1 || throw(ArgumentError("surface_PAR_division does not sum to 1"))
    
    fields = [CenterField(grid; 
                boundary_conditions = 
                    regularize_field_boundary_conditions(
                        FieldBoundaryConditions(top = ValueBoundaryCondition(ScaledSurfaceFunction(surface_PAR, surface_PAR_division[n]))), grid, name)) 
              for (n, name) in enumerate(field_names)]

    total_PAR = sum(fields) # need todo this before we convert `fields` to CuArray

    arch = architecture(grid)

    return MultiBandPhotosyntheticallyActiveRadiation(total_PAR,
                                                      on_architecture(arch, fields), 
                                                      field_names, 
                                                      on_architecture(arch, kʷ), 
                                                      on_architecture(arch, e), 
                                                      on_architecture(arch, χ), 
                                                      surface_PAR, 
                                                      on_architecture(arch, surface_PAR_division))
end

function numerical_mean(λ, C, idx1, idx2)
    ∫Cdλ = sum([(C[n] + C[n-1]) * (λ[n] - λ[n-1]) / 2 for n = idx1+1:idx2])
    ∫dλ = λ[idx2] - λ[idx1]
    return ∫Cdλ / ∫dλ
end

function par_symbol(n)
    subscripts = Symbol[]
    for digit in reverse(digits(n))
        push!(subscripts, Symbol(Char('\xe2\x82\x80'+digit)))
    end

    return Symbol(:PAR, subscripts...)
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
        @inbounds fields[n][i, j, k] = @inbounds surface_PAR * PAR_model.surface_PAR_division[n] * exp(zᶜ[grid.Nz] * (kʷ[n] + χ[n] * Chl[i, j, k] ^ e[n]))
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
    merge(NamedTuple{tuple(par.field_names...)}(par.fields), (PAR = par.total, ))

Adapt.adapt_structure(to, par::MultiBandPhotosyntheticallyActiveRadiation) = 
    MultiBandPhotosyntheticallyActiveRadiation(adapt(to, par.fields),
                                               nothing, # don't need this in the kernel
                                               adapt(to, par.water_attenuation_coefficient),
                                               adapt(to, par.chlorophyll_exponent),
                                               adapt(to, par.chlorophyll_attenuation_coefficient),
                                               adapt(to, par.surface_PAR),
                                               adapt(to, par.surface_PAR_division))
