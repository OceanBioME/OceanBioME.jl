struct TwoBandPhotosyntheticallyActiveRadiation{FT, FF, F, SPAR} <: AbstractSingleBandExponentialLightAttenuation{2, FF, F, SPAR}
    water_red_attenuation :: FT
    water_blue_attenuation :: FT
    chlorophyll_red_attenuation :: FT
    chlorophyll_blue_attenuation :: FT
    chlorophyll_red_exponent :: FT
    chlorophyll_blue_exponent :: FT
    pigment_ratio :: FT

    field :: F
    interface_field :: FF

    surface_PAR :: SPAR

    TwoBandPhotosyntheticallyActiveRadiation(water_red_attenuation::FT,
                                             water_blue_attenuation::FT,
                                             chlorophyll_red_attenuation::FT,
                                             chlorophyll_blue_attenuation::FT,
                                             chlorophyll_red_exponent::FT,
                                             chlorophyll_blue_exponent::FT,
                                             pigment_ratio::FT,
                                             field::F,
                                             interface_field::FF,
                                             surface_PAR::SPAR) where {FT, F, FF, SPAR} =
        new{FT, FF, F, SPAR}(water_red_attenuation,
                         water_blue_attenuation,
                         chlorophyll_red_attenuation,
                         chlorophyll_blue_attenuation,
                         chlorophyll_red_exponent,
                         chlorophyll_blue_exponent,
                         pigment_ratio,
                         field,
                         interface_field,
                         surface_PAR)
end

const TwoBandLight{FT, F, SPAR} = TwoBandPhotosyntheticallyActiveRadiation{FT, F, SPAR}

@inline attenuation(i, j, k, grid, la::TwoBandLight, clock, Chl, ::Val{1}) =
    @inbounds la.water_red_attenuation + (Chl[i, j, k]/la.pigment_ratio)^la.chlorophyll_red_exponent * la.chlorophyll_red_attenuation

@inline attenuation(i, j, k, grid, la::TwoBandLight, clock, Chl, ::Val{2}) =
    @inbounds la.water_blue_attenuation + (Chl[i, j, k]/la.pigment_ratio)^la.chlorophyll_blue_exponent * la.chlorophyll_blue_attenuation

"""
    TwoBandPhotosyntheticallyActiveRadiation(; grid::AbstractGrid{FT},
                                               water_red_attenuation = 0.225, # 1/m
                                               water_blue_attenuation = 0.0232, # 1/m
                                               chlorophyll_red_attenuation = 0.037, # 1/(m * (mgChl/m³) ^ eʳ)
                                               chlorophyll_blue_attenuation = 0.074, # 1/(m * (mgChl/m³) ^ eᵇ)
                                               chlorophyll_red_exponent = 0.629,
                                               chlorophyll_blue_exponent = 0.674,
                                               pigment_ratio = 0.7,
                                               surface_PAR = default_surface_PAR)

Keyword Arguments
==================

- `grid`: grid for building the model on
- `water_red_attenuation`, ..., `pigment_ratio`: parameter values
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
function TwoBandPhotosyntheticallyActiveRadiation(grid::AbstractGrid{FT}, surface_PAR;
                                                  water_red_attenuation = 0.225, # 1/m
                                                  water_blue_attenuation = 0.0232, # 1/m
                                                  chlorophyll_red_attenuation = 0.037, # 1/(m * (mgChl/m³) ^ eʳ)
                                                  chlorophyll_blue_attenuation = 0.074, # 1/(m * (mgChl/m³) ^ eᵇ)
                                                  chlorophyll_red_exponent = 0.629,
                                                  chlorophyll_blue_exponent = 0.674,
                                                  pigment_ratio = 0.7,
                                                  discrete_form = false,
                                                  parameters = nothing,
                                                  interface_field = nothing) where FT

    water_red_attenuation = convert(FT, water_red_attenuation)
    water_blue_attenuation = convert(FT, water_blue_attenuation)
    chlorophyll_red_attenuation = convert(FT, chlorophyll_red_attenuation)
    chlorophyll_blue_attenuation = convert(FT, chlorophyll_blue_attenuation)
    chlorophyll_red_exponent = convert(FT, chlorophyll_red_exponent)
    chlorophyll_blue_exponent = convert(FT, chlorophyll_blue_exponent)
    pigment_ratio = convert(FT, pigment_ratio)

    boundary_condition_kwargs = surface_PAR isa Function ? (; parameters, discrete_form) : NamedTuple()

    field = CenterField(grid; boundary_conditions =
                            regularize_field_boundary_conditions(
                                FieldBoundaryConditions(top = ValueBoundaryCondition(surface_PAR; boundary_condition_kwargs...)), grid, :PAR))

    # wrap surface_PAR to make it work with the `getbc` interface
    surface_PAR = materialize_condition(surface_PAR, parameters, discrete_form, ()) 
    surface_PAR = regularize_boundary_condition(surface_PAR, grid, (Center(), Center(), Center()), 3, RightBoundary, nothing)

    return TwoBandPhotosyntheticallyActiveRadiation(water_red_attenuation,
                                                    water_blue_attenuation,
                                                    chlorophyll_red_attenuation,
                                                    chlorophyll_blue_attenuation,
                                                    chlorophyll_red_exponent,
                                                    chlorophyll_blue_exponent,
                                                    pigment_ratio,
                                                    field,
                                                    interface_field,
                                                    surface_PAR)
end

summary(::TwoBandPhotosyntheticallyActiveRadiation{FT}) where {FT} = string("TwoBandPhotosyntheticallyActiveRadiation{$FT}")
show(io::IO, model::TwoBandPhotosyntheticallyActiveRadiation{FT}) where {FT} = print(io, summary(model))

biogeochemical_auxiliary_fields(par::TwoBandPhotosyntheticallyActiveRadiation) = (PAR = par.field, PAR_interface = par.interface_field)
biogeochemical_auxiliary_fields(par::TwoBandPhotosyntheticallyActiveRadiation{<:Any, Nothing}) = (PAR = par.field, )

adapt_structure(to, par::TwoBandPhotosyntheticallyActiveRadiation) =
    TwoBandPhotosyntheticallyActiveRadiation(adapt(to, par.water_red_attenuation),
                                             adapt(to, par.water_blue_attenuation),
                                             adapt(to, par.chlorophyll_red_attenuation),
                                             adapt(to, par.chlorophyll_blue_attenuation),
                                             adapt(to, par.chlorophyll_red_exponent),
                                             adapt(to, par.chlorophyll_blue_exponent),
                                             adapt(to, par.pigment_ratio),
                                             adapt(to, par.field),
                                             adapt(to, par.surface_PAR))
