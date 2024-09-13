function update_biogeochemical_state!(model, bgc::PISCES)
    # this should come from utils
    #update_mixed_layer_depth!(bgc, model)

    PAR = biogeochemical_auxiliary_fields(model.biogeochemistry.light_attenuation).PAR

    compute_euphotic_depth!(bgc.euphotic_depth, PAR)

    compute_mean_mixed_layer_vertical_diffusivity!(bgc.mean_mixed_layer_vertical_diffusivity, bgc.mixed_layer_depth, model)

    compute_mean_mixed_layer_light!(bgc.mean_mixed_layer_light, bgc.mixed_layer_depth, PAR, model)
    
    compute_calcite_saturation!(bgc.carbon_chemistry, bgc.calcite_saturation, model)

    #update_yearly_maximum_silicate!(bgc, model)

    return nothing
end