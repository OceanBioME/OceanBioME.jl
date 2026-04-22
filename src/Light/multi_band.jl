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
    MultiBandPhotosyntheticallyActiveRadiation(; grid::AbstractGrid{FT}, 
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
- `surface_PAR`: the photosynthetically available radiation at the surface, by default,
   assumed to have the 'continuous form' `condition(x, y t)`

If `parameters` is not `nothing`, then `surface_PAR` has the form
`func(x, y, t, parameters)`.

If `discrete_form = true`, surface_PAR is assumed to have the "discrete form",
```
condition(i, j, grid, clock, fields)
```
where `i`, and `j` are indices that vary along the boundary. If `discrete_form = true` and
`parameters` is not `nothing`, the function `condition` is called with
```
condition(i, j, grid, clock, fields, parameters)
```


"""
function MultiBandPhotosyntheticallyActiveRadiation(; grid::AbstractGrid{FT}, 
                                                      bands = ((400, 500), (500, 600), (600, 700)), #nm
                                                      base_bands = MOREL_λ,
                                                      base_water_attenuation_coefficient = MOREL_kʷ,
                                                      base_chlorophyll_exponent = MOREL_e,
                                                      base_chlorophyll_attenuation_coefficient = MOREL_χ,
                                                      field_names = ntuple(n->par_symbol(n), Val(length(bands))),
                                                      surface_PAR = default_surface_PAR,
                                                      discrete_form = false,
                                                      parameters = nothing,
                                                      surface_PAR_division = fill(1 / length(bands), length(bands)),) where FT
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
    
    n_fields = length(field_names)

    surface_boundary_conditions = 
        ntuple(n-> ValueBoundaryCondition(ScaledSurfaceFunction(surface_PAR, surface_PAR_division[n])), Val(n_fields))

    field_boundary_conditions =
        ntuple(n -> regularize_field_boundary_conditions(FieldBoundaryConditions(top = surface_boundary_conditions[n]),
                                                         grid, field_names[n]),
               Val(n_fields))

    fields = NamedTuple{field_names}(ntuple(n -> CenterField(grid; boundary_conditions = field_boundary_conditions[n]), Val(n_fields)))

    total_PAR = sum(fields)

    surface_PAR = materialize_condition(surface_PAR, parameters, discrete_form, ()) 
    surface_PAR = regularize_boundary_condition(surface_PAR, grid, (Center(), Center(), Center()), 3, RightBoundary, nothing)

    return MultiBandPhotosyntheticallyActiveRadiation(total_PAR,
                                                      fields, 
                                                      tuple(field_names...), 
                                                      convert.(FT, kʷ), 
                                                      convert.(FT, e), 
                                                      convert.(FT, χ), 
                                                      surface_PAR, 
                                                      surface_PAR_division)
end

function numerical_mean(λ, C, idx1, idx2)
    ∫Cdλ = sum([(C[n] + C[n-1]) * (λ[n] - λ[n-1]) / 2 for n = idx1+1:idx2])
    ∫dλ = λ[idx2] - λ[idx1]
    return ∫Cdλ / ∫dλ
end

@inline par_symbol(n) = Symbol(:PAR, Char('\xe2\x82\x80'+n)) #Symbol(:PAR, number_subscript(tuple(reverse(digits(n))...))...)

@inline number_subscript(digits::NTuple{N}) where N =
    ntuple(n->Symbol(Char('\xe2\x82\x80'+digits[n])), Val(N))

@kernel function update_MultiBandPhotosyntheticallyActiveRadiation!(grid, clock, field, kʷ, e, χ,
                                                                    _surface_PAR, 
                                                                    Chl, k′) 
    i, j = @index(Global, NTuple)

    surface_PAR = getbc(_surface_PAR, i, j, grid, clock, field)

    zᶜ = znodes(grid, Center(), Center(), Center())

    # first point below surface
    k = grid.Nz
    @inbounds field[i, j, k] = surface_PAR * exp(zᶜ[grid.Nz] * (kʷ + χ * Chl[i, j, k] ^ e))

    # the rest of the points
    for k in grid.Nz-1:-1:1
        Δz = @inbounds zᶜ[k] - zᶜ[k + 1] 
        @inbounds field[i, j, k] = @inbounds field[i, j, k + 1] * exp(Δz * (kʷ + χ * Chl[i, j, k] ^ e))
    end
end

function update_biogeochemical_state!(model, PAR::MultiBandPhotosyntheticallyActiveRadiation)
    grid = model.grid

    arch = architecture(grid)

    _, k′ = domain_boundary_indices(RightBoundary(), grid.Nz)

    for (n, field) in enumerate(PAR.fields)
        launch!(arch, grid, :xy, 
                update_MultiBandPhotosyntheticallyActiveRadiation!, grid, model.clock,
                field,
                PAR.water_attenuation_coefficient[n], 
                PAR.chlorophyll_exponent[n], 
                PAR.chlorophyll_attenuation_coefficient[n], 
                PAR.surface_PAR, 
                PAR.surface_PAR_division[n], 
                chlorophyll(model.biogeochemistry, model), 
                k′)
    end

    #for field in PAR.fields
    #    fill_halo_regions!(field, model.clock, fields(model))
    #end
end

summary(par::MultiBandPhotosyntheticallyActiveRadiation) = 
    string("Multi band light attenuation model with $(length(par.fields)) bands $(par.field_names)")

show(io::IO, model::MultiBandPhotosyntheticallyActiveRadiation) = print(io, summary(model))

biogeochemical_auxiliary_fields(par::MultiBandPhotosyntheticallyActiveRadiation) = 
    merge((PAR = par.total, ), par.fields)#NamedTuple{field_names(par.field_names, par.fields)}(par.fields))

@inline field_names(field_names, fields) = field_names
@inline field_names(::Nothing, fields::NTuple{N}) where N = ntuple(n -> par_symbol(n), Val(N))

# avoid passing this into kernels
Adapt.adapt_structure(to, par::MultiBandPhotosyntheticallyActiveRadiation) = 
    MultiBandPhotosyntheticallyActiveRadiation(adapt(to, par.total),
                                               adapt(to, par.fields),
                                               nothing,
                                               nothing, 
                                               nothing,
                                               nothing,
                                               nothing,
                                               nothing)

