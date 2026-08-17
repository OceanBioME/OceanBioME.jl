@kernel function update_TwoBandPhotosyntheticallyActiveRadiation!(PAR, grid, clock, Chl, surface_PAR, PAR_model)
    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)

    kʳ = PAR_model.water_red_attenuation
    kᵇ = PAR_model.water_blue_attenuation
    χʳ = PAR_model.chlorophyll_red_attenuation
    χᵇ = PAR_model.chlorophyll_blue_attenuation
    eʳ = PAR_model.chlorophyll_red_exponent
    eᵇ = PAR_model.chlorophyll_blue_exponent
    r  = PAR_model.pigment_ratio

    zᶜ = znodes(grid, Center(), Center(), Center())
    zᶠ = znodes(grid, Center(), Center(), Face())

    # first point below surface
    @inbounds begin
        ∫chlʳ = (zᶠ[grid.Nz + 1] - zᶜ[grid.Nz]) * (Chl[i, j, grid.Nz] / r)^eʳ
        ∫chlᵇ = (zᶠ[grid.Nz + 1] - zᶜ[grid.Nz]) * (Chl[i, j, grid.Nz] / r)^eᵇ
        PAR[i, j, grid.Nz] =  PAR⁰ * (exp(kʳ * zᶜ[grid.Nz] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[grid.Nz] - χᵇ * ∫chlᵇ)) / 2
    end

    # the rest of the points
    for k in grid.Nz-1:-1:1
        @inbounds begin
            ∫chlʳ += (zᶜ[k + 1] - zᶠ[k + 1]) * (Chl[i, j, k + 1] / r)^eʳ + (zᶠ[k + 1] - zᶜ[k]) * (Chl[i, j, k] / r)^eʳ
            ∫chlᵇ += (zᶜ[k + 1] - zᶠ[k + 1]) * (Chl[i, j, k + 1] / r)^eᵇ + (zᶠ[k + 1] - zᶜ[k]) * (Chl[i, j, k] / r)^eᵇ
            PAR[i, j, k] =  PAR⁰ * (exp(kʳ * zᶜ[k] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[k] - χᵇ * ∫chlᵇ)) / 2
        end
    end
end

struct TwoBandPhotosyntheticallyActiveRadiation{FT, F, SPAR}
    water_red_attenuation :: FT
    water_blue_attenuation :: FT
    chlorophyll_red_attenuation :: FT
    chlorophyll_blue_attenuation :: FT
    chlorophyll_red_exponent :: FT
    chlorophyll_blue_exponent :: FT
    pigment_ratio :: FT

    field :: F

    surface_PAR :: SPAR

    TwoBandPhotosyntheticallyActiveRadiation(water_red_attenuation::FT,
                                             water_blue_attenuation::FT,
                                             chlorophyll_red_attenuation::FT,
                                             chlorophyll_blue_attenuation::FT,
                                             chlorophyll_red_exponent::FT,
                                             chlorophyll_blue_exponent::FT,
                                             pigment_ratio::FT,
                                             field::F,
                                             surface_PAR::SPAR) where {FT, F, SPAR} =
        new{FT, F, SPAR}(water_red_attenuation,
                         water_blue_attenuation,
                         chlorophyll_red_attenuation,
                         chlorophyll_blue_attenuation,
                         chlorophyll_red_exponent,
                         chlorophyll_blue_exponent,
                         pigment_ratio,
                         field,
                         surface_PAR)
end

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

Total chlorophyll is obtained from `chlorophyll(model.biogeochemistry, model)`, allowing the
biogeochemistry to define chlorophyll from one or more plankton or chlorophyll fields.

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
                                                  parameters = nothing) where FT

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
                                                    surface_PAR)
end

function update_biogeochemical_state!(model, PAR::TwoBandPhotosyntheticallyActiveRadiation)
    arch = architecture(model.grid)
    Chl = chlorophyll(model.biogeochemistry, model)

    launch!(arch, model.grid, :xy, update_TwoBandPhotosyntheticallyActiveRadiation!, PAR.field, model.grid, model.clock, Chl, PAR.surface_PAR, PAR)
    #fill_halo_regions!(PAR.field, model.clock, fields(model))

    return nothing
end

summary(::TwoBandPhotosyntheticallyActiveRadiation{FT}) where {FT} = string("TwoBandPhotosyntheticallyActiveRadiation{$FT}")
show(io::IO, model::TwoBandPhotosyntheticallyActiveRadiation{FT}) where {FT} = print(io, summary(model))

biogeochemical_auxiliary_fields(par::TwoBandPhotosyntheticallyActiveRadiation) = (PAR = par.field, )

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
