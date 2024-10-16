# this file sets up the default configuration of Z and M which graze on P, D, (Z, ) and POC
function MicroAndMesoZooplankton(FT = Float64; 
    micro::MU = QualityDependantZooplankton(FT; maximum_grazing_rate = 3/day,
                                            food_preferences = (P = 1.0, D = 0.5, POC = 0.1, Z = 0),
                                            quadratic_mortality = 0.004/day,
                                            linear_mortality = 0.03/day,
                                            minimum_growth_efficiency = 0.3,
                                            maximum_flux_feeding_rate = 0.0,
                                            undissolved_calcite_fraction = 0.5,
                                            iron_ratio = 0.01),
    meso::ME = QualityDependantZooplankton(FT; maximum_grazing_rate = 0.75/day,
                                           food_preferences = (P = 0.3, D = 1.0, POC = 0.3, Z = 1.0),
                                           quadratic_mortality = 0.03/day,
                                           linear_mortality = 0.005/day,
                                           minimum_growth_efficiency = 0.35,
                                           maximum_flux_feeding_rate = 2e3 / 1e6, 
                                           undissolved_calcite_fraction = 0.75,
                                           iron_ratio = 0.015)) where {MU, ME}

    return MicroAndMeso{MU, ME, FT}(; micro, meso)
end

@inline concentration(::Val{:P},   i, j, k, fields) = @inbounds   fields.P[i, j, k]
@inline concentration(::Val{:D},   i, j, k, fields) = @inbounds   fields.D[i, j, k]
@inline concentration(::Val{:Z},   i, j, k, fields) = @inbounds   fields.Z[i, j, k]
@inline concentration(::Val{:POC}, i, j, k, fields) = @inbounds fields.POC[i, j, k]

@inline iron_ratio(::Val{:P},   i, j, k, bgc, fields) = @inbounds fields.PFe[i, j, k] / (fields.P[i, j, k] + eps(0.0))
@inline iron_ratio(::Val{:D},   i, j, k, bgc, fields) = @inbounds fields.DFe[i, j, k] / (fields.D[i, j, k] + eps(0.0))
@inline iron_ratio(::Val{:Z},   i, j, k, bgc, fields) = @inbounds bgc.zooplankton.micro.iron_ratio
@inline iron_ratio(::Val{:POC}, i, j, k, bgc, fields) = @inbounds fields.SFe[i, j, k] / (fields.POC[i, j, k] + eps(0.0))

@inline grazing_preference(val_prey_name, preferences) = 0
@inline grazing_preference(::Val{:P},     preferences) = preferences.P
@inline grazing_preference(::Val{:D},     preferences) = preferences.D
@inline grazing_preference(::Val{:Z},     preferences) = preferences.Z
@inline grazing_preference(::Val{:POC},   preferences) = preferences.POC

# TODO: this should dispatch on PISCES{<:NanoAndDiatoms, <:MicroAndMeso, <:Any, <:Two...} but the phyto and POC
# classes are not yet defined
@inline prey_names(::PISCES{<:Any, <:MicroAndMeso}, ::Val{:Z}) = (:P, :D, :POC)
@inline prey_names(::PISCES{<:Any, <:MicroAndMeso}, ::Val{:M}) = (:P, :D, :Z, :POC)

# TODO: move these somewhere else so they can be dispatched on ::PISCES{<:NanoAndDiatoms, <:MicroAndMeso, <:Any, <:TwoCompartementCarbonIronParticles}
@inline function extract_food_availability(::PISCES, i, j, k, fields, ::NTuple{N}) where N
    P = @inbounds fields.P[i, j, k]
    D = @inbounds fields.D[i, j, k]
    POC = @inbounds fields.POC[i, j, k]
    Z = @inbounds fields.Z[i, j, k]

    return (; P, D, POC, Z)
end

@inline function extract_iron_availability(bgc::PISCES, i, j, k, fields, ::NTuple{N}) where N
    P = @inbounds fields.PFe[i, j, k] / (fields.P[i, j, k] + eps(0.0))
    D = @inbounds fields.DFe[i, j, k] / (fields.D[i, j, k] + eps(0.0))
    POC = @inbounds fields.SFe[i, j, k] / (fields.POC[i, j, k] + eps(0.0))
    Z = bgc.zooplankton.micro.iron_ratio

    return (; P, D, POC, Z)
end