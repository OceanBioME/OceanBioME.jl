"""
    PhytoZoo

`PhytoZoo` defines the default living components for the `BiologyNutrientsDetritus` biogeochemical model.
It includes single `P`hytoplankton and `Z`ooplankton tracers which track the 
nitrogen (mmol N / m³). 

The phytoplankton evolve like `∂ₜP = (1-γ)μP - Gₚ - mP²` where μ is the growth rate (1/s)
which is of the form `phytoplankton_maximum_growth_rate` × light limitation × nutrient 
limitation. The `nutrient_limitation` term can be extended to account for different 
`nutrient` formulations. The `(1-γ)` term is the assimilated fraction. `Gₚ` accounts for 
grazing by the zooplankton, and `m` is the quadratic mortality rate (1/s/mmol N/m³).

The zooplankton evolves like `∂ₜZ = αG - mZ² - μZ` where `G` is the total grazing of 
phytoplankton and detritus with `α` the fraction assimilated, `m` is the quadratic 
mortality rate, and `μ` is the linear mortality rate taken to represent the excretion 
rate.
""" 
@kwdef struct PhytoZoo{FT}
            nitrate_half_saturation :: FT = 0.7       # mmol N/m³
            ammonia_half_saturation :: FT = 0.001     # mmol N/m³
               iron_half_saturation :: FT = 2e-4      # mmol Fe/m³ - estimated from fitting OSP 
         nitrate_ammonia_inhibition :: FT = 3.0       #
              light_half_saturation :: FT = 33.0      # W/m² 
  phytoplankton_maximum_growth_rate :: FT = 2.42e-5   # 1/s (double documented due to change to nitrogen limitation)
                         iron_ratio :: FT = 4.6375e-5 # mol Fe / mol N - from PISCES optimal ratio
   phytoplankton_exudation_fraction :: FT = 0.05      #
        ammonia_fraction_of_exudate :: FT = 0.75      #

            temperature_coefficient :: FT = 1.0       # Q10 factor, off by default 1.88 option from Kuhn, 2015

       phytoplankton_mortality_rate :: FT = 5.8e-7    # 1/s/mmol N/m³
         zooplankton_mortality_rate :: FT = 2.31e-6   # 1/s/mmol N/m³
         zooplankton_excretion_rate :: FT = 5.8e-7    # 1/s

       excretion_inorganic_fraction :: FT = 0.5       #
       preference_for_phytoplankton :: FT = 0.5       # 
               maximum_grazing_rate :: FT = 9.26e-6   # 1/s
            grazing_half_saturation :: FT = 1.0       # mmol N/m²
  zooplankton_assimilation_fraction :: FT = 0.7

    zooplankton_calcite_dissolution :: FT = 0.3

                     redfield_ratio :: FT = 6.56      # mol C/mol N
               carbon_calcate_ratio :: FT = 0.1       # mol CaCO₃/mol C
zooplankton_gut_calcite_dissolution :: FT = 0.3

    phytoplankton_chlorophyll_ratio :: FT = 1.31      # g Chl/mol N
end

const PHYTO_ZOO_BND = BiologyNutrientsDetritus{<:PhytoZoo}

required_biogeochemical_tracers(::PhytoZoo) = (:P, :Z)

@inline function (bnd::PHYTO_ZOO_BND)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    γ = bnd.biology.phytoplankton_exudation_fraction
    m = bnd.biology.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    μP = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    Gp = grazing(bnd, i, j, k, val_name, fields, auxiliary_fields)

    return (1 - γ) * μP - Gp - m * P^2
end # done!

@inline function (bnd::PHYTO_ZOO_BND)(i, j, k, grid, ::Val{:Z}, clock, fields, auxiliary_fields)
    α = bnd.biology.zooplankton_assimilation_fraction
    m = bnd.biology.zooplankton_mortality_rate
    μ = bnd.biology.zooplankton_excretion_rate

    Z = @inbounds fields.Z[i, j, k]
    G = total_grazing(bnd, i, j, k, fields, auxiliary_fields)

    return α * G - m * Z^2 - μ * Z
end # done!

@inline function total_grazing(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kG = bnd.biology.grazing_half_saturation
    g = bnd.biology.maximum_grazing_rate

    @inbounds begin
        Z = @inbounds fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(bnd.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(bnd.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    return g * food_concentration / (kG + food_concentration) * Z
end

@inline grazing_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) = 
    (1 - bnd.biology.zooplankton_assimilation_fraction) * total_grazing(bnd, i, j, k, fields, auxiliary_fields)

@inline function phytoplankton_preference(bnd, i, j, k, fields, auxiliary_fields)
    p̃ = bnd.biology.preference_for_phytoplankton

    P = @inbounds fields.P[i, j, k]
    sPOM = small_particulate_concentration(bnd.detritus, i, j, k, fields, auxiliary_fields)
    
    return p̃ * P / (p̃ * P + (1 - p̃) * sPOM + eps(0.0))
end

##### growth limitation
@inline function nitrogen_limitation(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kNO₃ = bnd.biology.nitrate_half_saturation
    kNH₄ = bnd.biology.ammonia_half_saturation
    ψ = bnd.biology.nitrate_ammonia_inhibition

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    # factor of 1/2 so that the limitation varies between 0 and 1
    # correspondingly μ₀ is doubled
    return (nitrate_limitation + ammonia_limitation) / 2
end

@inline nutrient_limitation(bnd::LOBSTER{<:Any, <:NitrateAmmonia}, i, j, k, fields, auxiliary_fields) =
    nitrogen_limitation(bnd, i, j, k, fields, auxiliary_fields)

@inline function nutrient_limitation(bnd::LOBSTER{<:Any, <:NitrateAmmoniaIron}, i, j, k, fields, auxiliary_fields)
    kFe = bnd.biology.iron_half_saturation
    Fe = @inbounds fields.Fe[i, j, k]

    iron_limitation = Fe / (kFe + Fe)

    return nitrogen_limitation(bnd, i, j, k, fields, auxiliary_fields) * iron_limitation
end

@inline function nutrient_limitation(bnd::LOBSTER{<:Any, <:Nutrient}, i, j, k, fields, auxiliary_fields)
    kNO₃ = bnd.biology.nitrate_half_saturation

    @inbounds begin
        N = fields.N[i, j, k]
    end

    return N / (N + kNO₃)
end

@inline function nitrogen_limitation(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kNO₃ = bnd.biology.nitrate_half_saturation
    kNH₄ = bnd.biology.ammonia_half_saturation
    ψ = bnd.biology.nitrate_ammonia_inhibition

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    # factor of 1/2 so that the limitation varies between 0 and 1
    # correspondingly μ₀ is doubled
    return (nitrate_limitation + ammonia_limitation) / 2
end


##### total phytoplankton growth
@inline function phytoplankton_growth(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kPAR = bnd.biology.light_half_saturation
    μ₀   = bnd.biology.phytoplankton_maximum_growth_rate
    Q10  = bnd.biology.temperature_coefficient

    @inbounds begin
        PAR = auxiliary_fields.PAR[i, j, k]
        P   = fields.P[i, j, k]
        T   = possibly_get_T(i, j, k, fields)
    end

    Lₙ = nutrient_limitation(bnd, i, j, k, fields, auxiliary_fields)
    Lₗ = PAR / (kPAR + PAR)
    Lₜ = Q10 ^ (T / 10)

    return μ₀ * Lₗ * Lₙ * Lₜ * P
end

@inline function possibly_get_T(i, j, k, fields::NamedTuple{tracer_names}) where tracer_names
    if :T ∈ tracer_names # this is a compile time evaluation
        return fields.T[i, j, k]
    else
        return zero(eltype(fields.P)) # this doesn't seem optimal but we don't pass grid...
    end
end

###### nutrient uptake
@inline function nutrient_uptake(bnd::PHYTO_ZOO_BND, i, j, k, val_name::Val{:NO₃}, fields, auxiliary_fields)
    μ = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    kNO₃ = bnd.biology.nitrate_half_saturation
    kNH₄ = bnd.biology.ammonia_half_saturation
    ψ = bnd.biology.nitrate_ammonia_inhibition

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    return μ * nitrate_limitation / (nitrate_limitation + ammonia_limitation + eps(0.0))
end

@inline function nutrient_uptake(bnd::PHYTO_ZOO_BND, i, j, k, val_name::Val{:NH₄}, fields, auxiliary_fields)
    kNO₃ = bnd.biology.nitrate_half_saturation
    kNH₄ = bnd.biology.ammonia_half_saturation
    ψ    = bnd.biology.nitrate_ammonia_inhibition
    α    = bnd.biology.ammonia_fraction_of_exudate
    γ    = bnd.biology.phytoplankton_exudation_fraction
    
    μ = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    waste = α * γ * μ

    return μ * ammonia_limitation / (nitrate_limitation + ammonia_limitation + eps(0.0)) - waste
end

@inline function nutrient_uptake(bnd::PHYTO_ZOO_BND, i, j, k, val_name::Val{:Fe}, fields, auxiliary_fields)
    R = bnd.biology.iron_ratio
    
    μ = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    return R * μ
end

@inline nutrient_uptake(bnd::PHYTO_ZOO_BND, i, j, k, val_name::Val{:N}, fields, auxiliary_fields) =
    phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

@inline function phytoplankton_primary_production(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α = bnd.biology.ammonia_fraction_of_exudate
    γ = bnd.biology.phytoplankton_exudation_fraction
    ρ = bnd.biology.carbon_calcate_ratio
    R = bnd.biology.redfield_ratio

    μP = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    return (1 + ρ * (1 - γ) - α * γ) * μP * R
end

##### mortality waste
@inline function biology_inorganic_nitrogen_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α = bnd.biology.excretion_inorganic_fraction
    μ = bnd.biology.zooplankton_excretion_rate

    Z = @inbounds fields.Z[i, j, k]

    return α * μ * Z
end

@inline biology_inorganic_carbon_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    bnd.biology.redfield_ratio * biology_inorganic_nitrogen_waste(bnd, i, j, k, fields, auxiliary_fields)


@inline function biology_organic_nitrogen_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    αZ = bnd.biology.excretion_inorganic_fraction
    μ  = bnd.biology.zooplankton_excretion_rate
    αP = bnd.biology.ammonia_fraction_of_exudate
    γ  = bnd.biology.phytoplankton_exudation_fraction

    Z = @inbounds fields.Z[i, j, k]

    μP = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    return (1 - αP) * γ * μP + (1 - αZ) * μ * Z
end

@inline biology_organic_carbon_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    bnd.biology.redfield_ratio * biology_organic_nitrogen_waste(bnd, i, j, k, fields, auxiliary_fields)

@inline function solid_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α  = bnd.biology.zooplankton_assimilation_fraction
    mᴾ = bnd.biology.phytoplankton_mortality_rate
    mᶻ = bnd.biology.zooplankton_mortality_rate

    @inbounds begin
        P = fields.P[i, j, k]
        Z = fields.Z[i, j, k]
    end

    G = total_grazing(bnd, i, j, k, fields, auxiliary_fields)

    return (1 - α) * G + mᴾ * P^2 + mᶻ * Z^2
end

@inline solid_carbon_waste(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    solid_waste(bnd, i, j, k, fields, auxiliary_fields) * bnd.biology.redfield_ratio

#### grazing
@inline function grazing(bnd::PHYTO_ZOO_BND, i, j, k, ::Val{:P}, fields, auxiliary_fields)
    kG = bnd.biology.grazing_half_saturation
    g  = bnd.biology.maximum_grazing_rate

    @inbounds begin
        Z = fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(bnd.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(bnd.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    return g * p * P / (kG + food_concentration) * Z
end

@inline function grazing(bnd::PHYTO_ZOO_BND, i, j, k, ::Union{Val{:sPOM}, Val{:sPON}}, fields, auxiliary_fields)
    kG = bnd.biology.grazing_half_saturation
    g  = bnd.biology.maximum_grazing_rate

    @inbounds begin
        Z = fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(bnd.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(bnd.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    return g * (1-p) * sPOM / (kG + food_concentration) * Z
end

@inline function grazing(bnd::PHYTO_ZOO_BND, i, j, k, ::Val{:D}, fields, auxiliary_fields)
    kG = bnd.biology.grazing_half_saturation
    g  = bnd.biology.maximum_grazing_rate

    @inbounds begin
        Z = fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(bnd.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(bnd.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    return g * (1-p) * sPOM / (kG + food_concentration) * Z
end

@inline grazing(bnd::PHYTO_ZOO_BND, i, j, k, ::Val{:sPOC}, fields, auxiliary_fields) =
    grazing(bnd, i, j, k, Val(:sPOM), fields, auxiliary_fields) * bnd.biology.redfield_ratio

@inline function weighted_phytoplankton_preference(biology::PhytoZoo, P, sPOM)
    p̃ = biology.preference_for_phytoplankton

    return p̃ * P / (p̃ * P + (1 - p̃) * sPOM + eps(0.0))
end

@inline function calcite_production(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α  = bnd.biology.zooplankton_assimilation_fraction
    mᴾ = bnd.biology.phytoplankton_mortality_rate
    R  = bnd.biology.redfield_ratio 
    ρ  = bnd.biology.carbon_calcate_ratio
    η  = bnd.biology.zooplankton_gut_calcite_dissolution

    P = @inbounds fields.P[i, j, k]

    G = grazing(bnd, i, j, k, Val(:P), fields, auxiliary_fields)

    return (G * (1 - η) + mᴾ * P^2) * ρ * R
end

@inline function calcite_dissolution(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    R  = bnd.biology.redfield_ratio 
    ρ  = bnd.biology.carbon_calcate_ratio
    η  = bnd.biology.zooplankton_gut_calcite_dissolution

    G = grazing(bnd, i, j, k, Val(:P), fields, auxiliary_fields)

    return G * η * ρ * R
end

# when we have variable redfield particles the calcite production goes to bPOC, but this
# isn't implicitly captured in the bPOM when we have a fixed redfield, so we have to 
# assume instant dissolution and put it back in DIC (or we could choose to lose it),
# maybe that would be better?
# but when we don't this can't be implicitly captured so we put it all back in the DIC compartement
@inline function calcite_dissolution(bnd::LOBSTER{<:PhytoZoo, <:Any, <:TwoParticleAndDissolved}, i, j, k, fields, auxiliary_fields)
    R  = bnd.biology.redfield_ratio 
    ρ  = bnd.biology.carbon_calcate_ratio
    mᴾ = bnd.biology.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    G = grazing(bnd, i, j, k, Val(:P), fields, auxiliary_fields)

    return (G + mᴾ * P^2) * ρ * R
end

@inline function calcite_uptake(bnd::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    R  = bnd.biology.redfield_ratio 
    ρ  = bnd.biology.carbon_calcate_ratio

    μP = phytoplankton_growth(bnd, i, j, k, fields, auxiliary_fields)

    return 2 * ρ * μP * R
end