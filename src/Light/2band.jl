@kernel function update_TwoBandPhotosyntheticallyActiveRatiation!(PAR, grid, P, surface_PAR, t, PAR_model) 
    i, j = @index(Global, NTuple)
    x, y = xnode(Center(), i, grid), ynode(Center(), j, grid)
    
    PAR⁰ = surface_PAR(x, y, t)

    kʳ = PAR_model.water_red_attenuation
    kᵇ = PAR_model.water_blue_attenuation
    χʳ = PAR_model.chlorophyll_red_attenuation
    χᵇ = PAR_model.chlorophyll_blue_attenuation
    eʳ = PAR_model.chlorophyll_red_exponent
    eᵇ = PAR_model.chlorophyll_blue_exponent
    r = PAR_model.pigment_ratio
    Rᶜₚ = PAR_model.phytoplankton_chlorophyll_ratio

    zᶜ = znodes(Center, grid)
    zᶠ = znodes(Face, grid)
    
    ∫chlʳ = @inbounds (zᶠ[grid.Nz + 1] - zᶜ[grid.Nz]) * (P[i, j, grid.Nz] * Rᶜₚ / r) ^ eʳ
    ∫chlᵇ = @inbounds (zᶠ[grid.Nz + 1] - zᶜ[grid.Nz]) * (P[i, j, grid.Nz] * Rᶜₚ / r) ^ eᵇ

    # first point below surface
    @inbounds PAR[i, j, grid.Nz] =  PAR⁰ * (exp(kʳ * zᶜ[grid.Nz] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[grid.Nz] - χᵇ * ∫chlᵇ)) / 2

    # first point above surface (so that we interpolate to PAR⁰ at surface)
    @inbounds PAR[i, j, grid.Nz + 1] = PAR⁰ * (2 - (exp(kʳ * zᶜ[grid.Nz] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[grid.Nz] - χᵇ * ∫chlᵇ)) / 2)

    @inbounds for k in grid.Nz-1:-1:1
        ∫chlʳ += (zᶜ[k + 1] - zᶠ[k + 1]) * (P[i, j, k+1] * Rᶜₚ / r) ^ eʳ + (zᶠ[k + 1] - zᶜ[k]) * (P[i, j, k] * Rᶜₚ / r) ^ eʳ
        ∫chlᵇ += (zᶜ[k + 1] - zᶠ[k + 1]) * (P[i, j, k+1] * Rᶜₚ / r) ^ eᵇ + (zᶠ[k + 1] - zᶜ[k]) * (P[i, j, k] * Rᶜₚ / r) ^ eᵇ
        PAR[i, j, k] =  PAR⁰ * (exp(kʳ * zᶜ[k] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[k] - χᵇ * ∫chlᵇ)) / 2
    end
end 

struct TwoBandPhotosyntheticallyActiveRatiation{FT, F}
    water_red_attenuation :: FT
    water_blue_attenuation :: FT
    chlorophyll_red_attenuation :: FT
    chlorophyll_blue_attenuation :: FT
    chlorophyll_red_exponent :: FT
    chlorophyll_blue_exponent :: FT
    pigment_ratio :: FT

    phytoplankton_chlorophyll_ratio :: FT

    field :: F

    function TwoBandPhotosyntheticallyActiveRatiation(water_red_attenuation::FT,
                                                      water_blue_attenuation::FT,
                                                      chlorophyll_red_attenuation::FT,
                                                      chlorophyll_blue_attenuation::FT,
                                                      chlorophyll_red_exponent::FT,
                                                      chlorophyll_blue_exponent::FT,
                                                      pigment_ratio::FT,
                                                      phytoplankton_chlorophyll_ratio::FT,
                                                      field::F) where {FT, F}
        return new{FT, F}(water_red_attenuation,
                          water_blue_attenuation,
                          chlorophyll_red_attenuation,
                          chlorophyll_blue_attenuation,
                          chlorophyll_red_exponent,
                          chlorophyll_blue_exponent,
                          pigment_ratio,
                          phytoplankton_chlorophyll_ratio,
                          field)
    end
end

function TwoBandPhotosyntheticallyActiveRatiation(; grid, 
                                                    water_red_attenuation::FT = 0.225, # 1/m
                                                    water_blue_attenuation::FT = 0.0232, # 1/m
                                                    chlorophyll_red_attenuation::FT = 0.037, # 1/(m * (mgChl/m³) ^ eʳ)
                                                    chlorophyll_blue_attenuation::FT = 0.074, # 1/(m * (mgChl/m³) ^ eᵇ)
                                                    chlorophyll_red_exponent::FT = 0.629,
                                                    chlorophyll_blue_exponent::FT = 0.674,
                                                    pigment_ratio::FT = 0.7,
                                                    phytoplankton_chlorophyll_ratio::FT = 1.31) where FT # mgChl/mol N
    
    field = CenterField(grid; boundary_conditions = FieldBoundaryConditions(grid, (Center, Center, Center); 
                                                                            top = BoundaryCondition(PAR_boundary, nothing_function)))

    return TwoBandPhotosyntheticallyActiveRatiation(water_red_attenuation,
                                                    water_blue_attenuation,
                                                    chlorophyll_red_attenuation,
                                                    chlorophyll_blue_attenuation,
                                                    chlorophyll_red_exponent,
                                                    chlorophyll_blue_exponent,
                                                    pigment_ratio,
                                                    phytoplankton_chlorophyll_ratio,
                                                    field)
end


function update_PAR!(model, PAR::TwoBandPhotosyntheticallyActiveRatiation, surface_PAR)
    arch = architecture(model.grid)
    event = launch!(arch, model.grid, :xy, update_TwoBandPhotosyntheticallyActiveRatiation!, model.auxiliary_fields.PAR, model.grid, model.tracers.P, surface_PAR, model.clock.time, PAR)
    wait(event)
end

required_PAR_fields(::TwoBandPhotosyntheticallyActiveRatiation) = (:PAR, )

summary(::TwoBandPhotosyntheticallyActiveRatiation{FT}) where {FT} = string("Two-band light attenuation model ($FT)")
show(io::IO, model::TwoBandPhotosyntheticallyActiveRatiation{FT}) where {FT} = print(io, summary(model))

biogeochemical_auxiliary_fieilds(par::TwoBandPhotosyntheticallyActiveRatiation) = (PAR = par.field, )

adapt_structure(to, par::TwoBandPhotosyntheticallyActiveRatiation) = 
    TwoBandPhotosyntheticallyActiveRatiation(par.water_red_attenuation,
                                             par.water_blue_attenuation,
                                             par.chlorophyll_red_attenuation,
                                             par.chlorophyll_blue_attenuation,
                                             par.chlorophyll_red_exponent,
                                             par.chlorophyll_blue_exponent,
                                             par.pigment_ratio,
                                             par.phytoplankton_chlorophyll_ratio,
                                             adapt_structure(to, par.field))