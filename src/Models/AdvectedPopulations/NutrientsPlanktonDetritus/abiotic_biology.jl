
"""
    Abiotic

`Abiotic` defines a plankton with no tracers where production is
instantaneously transformed into detritus.

This is the simple model available in MITGCM [1] and described in [2]

[1]: https://mitgcm.readthedocs.io/en/latest/examples/global_oce_biogeo/global_oce_biogeo.html
[2]: https://web.mit.edu/globalchange/www/MITJPSPGC_Rpt122.pdf
"""

@kwdef struct Abiotic{FT}
       light_half_saturation :: FT = 25 # W/m²
    nutrient_half_saturation :: FT = 0.1 # mmolP/m³
          maximum_production :: FT = 0.3 / (365days) # mmolP/s - this seems very small
    dissolved_waste_fraction :: FT = 0.67
    carbon_to_nutrient_ratio :: FT = 117.0 # this is for phosphate
  nitrogen_to_nutrient_ratio :: FT = 16.0 # for alkalinity - can be 1 if the nutrient is nitrate
                  rain_ratio :: FT = 0.07
end

required_biogeochemical_tracers(::Abiotic) = () # not necessary but for clarity
required_biogeochemical_auxiliary_fields(::Abiotic) = (:PAR, ) # not necessary but for clarity

# MITGCM assumes that the solid fraction of production instantly sinks and is remineralised
# This means you don't have to have a particulate tracer, and you aren't limited by the 
# sinking CFL (which you are in coarse global runs, if you have 10m surface spacing and want a 
# time step of 1.5hours max speed is 80m/day)
# TODO: implement this, for now, I will also implement a simple dissolved particulate detritus

@inline function phytoplankton_primary_production(bgc::NutrientsPlanktonDetritus{<:Any, <:Abiotic{FT}}, i, j, k, fields, auxiliary_fields)
    μ    = bgc.plankton.maximum_production
    kPAR = bgc.plankton.light_half_saturation
    r    = bgc.plankton.carbon_to_nutrient_ratio

    @inbounds begin
        PAR = auxiliary_fields.PAR[i, j, k]
    end

    Ln = nutrient_limitation(bgc, i, j, k, fields, auxiliary_fields)
    Ll = PAR / (PAR + kPAR)

    return μ * Ln * Ll * r
end

@inline plankton_inorganic_carbon_waste(::NutrientsPlanktonDetritus{<:Any, <:Abiotic{FT}}, args...) where FT = zero(FT)

@inline nutrient_uptake(bgc::NutrientsPlanktonDetritus{<:Any, <:Abiotic{FT}}, args...) =
    phytoplankton_primary_production(bgc, args...) / bgc.plankton.carbon_to_nutrient_ratio

@inline plankton_organic_nitrogen_waste(bgc,::NutrientsPlanktonDetritus{<:Any, <:Abiotic}, args...) =
    bgc.plankton.dissolved_waste_fraction * phytoplankton_primary_production(bgc, args...) / bgc.plankton.carbon_to_nutrient_ratio

@inline solid_waste(bgc::NutrientsPlanktonDetritus{<:Any, <:Abiotic}, args...) =
    (1 - bgc.plankton.dissolved_waste_fraction) * phytoplankton_primary_production(bgc, args...) / bgc.plankton.carbon_to_nutrient_ratio

@inline grazing(::NutrientsPlanktonDetritus{<:Any, <:Abiotic{FT}}, args...) where FT = zero(FT)

@inline calcite_uptake(bgc::NutrientsPlanktonDetritus{<:Any, <:Abiotic}, args...) = 
    bgc.plankton.rain_ratio * (1 - bgc.plankton.dissolved_waste_fraction) * phytoplankton_primary_production(bgc, args...)

@inline calcite_dissolution(bgc::NutrientsPlanktonDetritus{<:Any, <:Abiotic}, agrs...) = 
    -calcite_uptake(bgc, args...)

# TODO: make this generic so we can have any number of limiting nutrients
@inline function nutrient_limitation(bgc::NutrientsPlanktonDetritus{<:Nutrient, <:Abiotic}, i, j, k, fields, auxiliary_fields)
    k = bgc.plankton.nutrient_half_saturation

    @inbounds begin
        N = fields.N[i, j, k]
    end

    return N / (N + k)
end
