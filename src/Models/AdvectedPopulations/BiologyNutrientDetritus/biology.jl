"""
    PhytoZoo

`PhytoZoo` defines the default living component for the `LOBSTER` biogeochemical model.
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
@kwdef struct PhytoZoo{FT, LL, TC, PM, GC, PS, ZS}
            nitrate_half_saturation :: FT = 0.7       # mmol N/m³
            ammonia_half_saturation :: FT = 0.001     # mmol N/m³
               iron_half_saturation :: FT = 2e-4      # mmol Fe/m³ - estimated from fitting OSP 
         nitrate_ammonia_inhibition :: FT = 3.0       #
              light_half_saturation :: FT = 33.0      # W/m² 
  phytoplankton_maximum_growth_rate :: FT = 2.42e-5   # 1/s (double documented due to change to nitrogen limitation)
                         iron_ratio :: FT = 4.6375e-5 # mol Fe / mol N - from PISCES optimal ratio
   phytoplankton_exudation_fraction :: FT = 0.05      #
        ammonia_fraction_of_exudate :: FT = 0.75      #
                   light_limitation :: LL = MondoLightLimitation()

            temperature_coefficient :: TC = nothing   # Q10 factor, off by default 1.88 option from Kuhn, 2015

       phytoplankton_mortality_rate :: FT = 5.8e-7    # 1/s/mmol N/m³
         zooplankton_mortality_rate :: FT = 2.31e-6   # 1/s/mmol N/m³
         zooplankton_excretion_rate :: FT = 5.8e-7    # 1/s

phytoplankton_mortality_formulation :: PM = Quadratic()

 phytoplankton_solid_waste_fraction :: FT = 1.0
       excretion_inorganic_fraction :: FT = 0.5       #
       preference_for_phytoplankton :: FT = 0.5       # 
               maximum_grazing_rate :: FT = 9.26e-6   # 1/s
            grazing_half_saturation :: FT = 1.0       # mmol N/m²
  zooplankton_assimilation_fraction :: FT = 0.7

  grazing_concentration_formulation :: GC = Quadratic()

    zooplankton_calcite_dissolution :: FT = 0.3

                     redfield_ratio :: FT = 6.56      # mol C/mol N
               carbon_calcate_ratio :: FT = 0.1       # mol CaCO₃/mol C
zooplankton_gut_calcite_dissolution :: FT = 0.3

    phytoplankton_chlorophyll_ratio :: FT = 1.31      # g Chl/mol N

     phytoplankton_sinking_velocity :: PS = (u = ZeroField(), v = ZeroField(), w = ZeroField()) # this is already the fallback in Oceananigans
       zooplankton_sinking_velocity :: ZS = (u = ZeroField(), v = ZeroField(), w = ZeroField())
end

biogeochemical_drift_velocity(phytozoo::PhytoZoo, ::Val{:P}) = 
    phytozoo.phytoplankton_sinking_velocity
biogeochemical_drift_velocity(phytozoo::PhytoZoo, ::Val{:Z}) = 
    phytozoo.zooplankton_sinking_velocity

function PhytoZoo(grid; 
                  phytoplankton_sinking_speed = 0.0, # m/s
                  zooplankton_sinking_speed = 0.0, # m/s
                  open_bottom = true,
                  kwargs...)


    phytoplankton_sinking_velocity = setup_velocity_fields((; P = phytoplankton_sinking_speed), grid, open_bottom; three_D = true).P
    zooplankton_sinking_velocity = setup_velocity_fields((; Z = zooplankton_sinking_speed), grid, open_bottom; three_D = true).Z

    return PhytoZoo(; phytoplankton_sinking_velocity, zooplankton_sinking_velocity, kwargs...)
end

const PHYTO_ZOO_BND = BiologyNutrientDetritus{<:Any, <:PhytoZoo}

required_biogeochemical_tracers(pz::PhytoZoo) = ((:P, :Z)..., (isnothing(pz.temperature_coefficient) ? () : (:T, ))...)
required_biogeochemical_auxiliary_fields(::PhytoZoo) = (:PAR, )

struct Linear end
struct Quadratic end

@inline mortality(::Linear, X, m) = m * X
@inline mortality(::Quadratic, X, m) = m * X^2

@inline concentration_limit(::Linear, X, k) = X / (X + k)
@inline concentration_limit(::Quadratic, X, k) = X^2 / (X^2 + k^2)

@inline function (lobster::PHYTO_ZOO_BND)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    γ = lobster.biology.phytoplankton_exudation_fraction
    m = lobster.biology.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    Gp = grazing(lobster, i, j, k, val_name, fields, auxiliary_fields)

    ν  = mortality(lobster.biology.phytoplankton_mortality_formulation, P, m)

    return (1 - γ) * μP - Gp - ν
end 

@inline function (lobster::PHYTO_ZOO_BND)(i, j, k, grid, ::Val{:Z}, clock, fields, auxiliary_fields)
    α = lobster.biology.zooplankton_assimilation_fraction
    m = lobster.biology.zooplankton_mortality_rate
    μ = lobster.biology.zooplankton_excretion_rate

    Z = @inbounds fields.Z[i, j, k]
    G = total_grazing(lobster, i, j, k, fields, auxiliary_fields)

    return α * G - m * Z^2 - μ * Z
end 

@inline function total_grazing(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kG = lobster.biology.grazing_half_saturation
    g = lobster.biology.maximum_grazing_rate

    @inbounds begin
        Z = @inbounds fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(lobster.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    L = concentration_limit(lobster.biology.grazing_concentration_formulation, food_concentration, kG)

    return g * L * Z
end

@inline grazing_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) = 
    (1 - lobster.biology.zooplankton_assimilation_fraction) * total_grazing(lobster, i, j, k, fields, auxiliary_fields)

@inline function phytoplankton_preference(lobster, i, j, k, fields, auxiliary_fields)
    p̃ = lobster.biology.preference_for_phytoplankton

    P = @inbounds fields.P[i, j, k]
    sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    
    return p̃ * P / (p̃ * P + (1 - p̃) * sPOM + eps(0.0))
end

##### growth limitation
@inline function nitrogen_limitation(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kNO₃ = lobster.biology.nitrate_half_saturation
    kNH₄ = lobster.biology.ammonia_half_saturation
    ψ = lobster.biology.nitrate_ammonia_inhibition

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

@inline nutrient_limitation(lobster::BiologyNutrientDetritus{<:NitrateAmmonia}, i, j, k, fields, auxiliary_fields) =
    nitrogen_limitation(lobster, i, j, k, fields, auxiliary_fields)

@inline function nutrient_limitation(lobster::BiologyNutrientDetritus{<:NitrateAmmoniaIron}, i, j, k, fields, auxiliary_fields)
    kFe = lobster.biology.iron_half_saturation
    Fe = @inbounds fields.Fe[i, j, k]

    iron_limitation = Fe / (kFe + Fe)

    return nitrogen_limitation(lobster, i, j, k, fields, auxiliary_fields) * iron_limitation
end


@inline function nutrient_limitation(lobster::BiologyNutrientDetritus{<:Nutrient}, i, j, k, fields, auxiliary_fields)
    kNO₃ = lobster.biology.nitrate_half_saturation

    @inbounds begin
        N = fields.N[i, j, k]
    end

    return N / (N + kNO₃)
end

##### total phytoplankton growth
@inline function phytoplankton_growth(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    kPAR = lobster.biology.light_half_saturation
    μ₀   = lobster.biology.phytoplankton_maximum_growth_rate
    Q10  = lobster.biology.temperature_coefficient

    @inbounds begin
        PAR = auxiliary_fields.PAR[i, j, k]
        P   = fields.P[i, j, k]
    end

    Lₙ = nutrient_limitation(lobster, i, j, k, fields, auxiliary_fields)
    Lₗ = light_limitation(lobster.biology.light_limitation, PAR, kPAR)
    Lₜ = temperature_limitation(Q10, i, j, k, fields)

    return μ₀ * Lₗ * Lₙ * Lₜ * P
end

@inline temperature_limitation(::Nothing, i, j, k, fields) = zero(eltype(fields.P))

@inline function temperature_limitation(Q10, i, j, k, fields)
    T = @inbounds fields.T[i, j, k]

    return Q10 ^ (T / 10)
end

struct MondoLightLimitation end
struct AnalyticalLightLimitation end # This formulation is justified because you can rearrange to αI/√(μ² + α²I²) which if you put simple I(t) in can be integrated analytically

@inline light_limitation(::MondoLightLimitation, PAR, kPAR) = PAR / (kPAR + PAR)
@inline light_limitation(::AnalyticalLightLimitation, PAR, kPAR) = PAR / sqrt(PAR^2 + kPAR^2)

###### nutrient uptake
@inline function nutrient_uptake(lobster::PHYTO_ZOO_BND, i, j, k, val_name::Val{:NO₃}, fields, auxiliary_fields)
    μ = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    kNO₃ = lobster.biology.nitrate_half_saturation
    kNH₄ = lobster.biology.ammonia_half_saturation
    ψ = lobster.biology.nitrate_ammonia_inhibition

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    return μ * nitrate_limitation / (nitrate_limitation + ammonia_limitation + eps(0.0))
end

@inline function nutrient_uptake(lobster::PHYTO_ZOO_BND, i, j, k, val_name::Val{:NH₄}, fields, auxiliary_fields)
    kNO₃ = lobster.biology.nitrate_half_saturation
    kNH₄ = lobster.biology.ammonia_half_saturation
    ψ    = lobster.biology.nitrate_ammonia_inhibition
    α    = lobster.biology.ammonia_fraction_of_exudate
    γ    = lobster.biology.phytoplankton_exudation_fraction
    
    μ = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    waste = α * γ * μ

    return μ * ammonia_limitation / (nitrate_limitation + ammonia_limitation + eps(0.0)) - waste
end

@inline function nutrient_uptake(lobster::PHYTO_ZOO_BND, i, j, k, val_name::Val{:Fe}, fields, auxiliary_fields)
    R = lobster.biology.iron_ratio
    
    μ = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    return R * μ
end

@inline function nutrient_uptake(lobster::PHYTO_ZOO_BND, i, j, k, val_name::Val{:N}, fields, auxiliary_fields)
    α = lobster.biology.ammonia_fraction_of_exudate
    γ = lobster.biology.phytoplankton_exudation_fraction
    
    μ = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    return μ * (1 - α * γ)
end

@inline function phytoplankton_primary_production(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α = lobster.biology.ammonia_fraction_of_exudate
    γ = lobster.biology.phytoplankton_exudation_fraction
    ρ = lobster.biology.carbon_calcate_ratio
    R = lobster.biology.redfield_ratio

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    return (1 + ρ * (1 - γ) - α * γ) * μP * R
end

##### mortality waste
@inline function biology_inorganic_nitrogen_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    αᴾ = lobster.biology.phytoplankton_solid_waste_fraction
    αᶻ = lobster.biology.excretion_inorganic_fraction
    μ  = lobster.biology.zooplankton_excretion_rate
    mᴾ = lobster.biology.phytoplankton_mortality_rate

    @inbounds begin
        P = fields.P[i, j, k]
        Z = fields.Z[i, j, k]
    end

    νP = mortality(lobster.biology.phytoplankton_mortality_formulation, P, mᴾ)

    return αᶻ * μ * Z + (1 - αᴾ) * νP
end

@inline biology_inorganic_carbon_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    lobster.biology.redfield_ratio * biology_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)


@inline function biology_organic_nitrogen_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    αZ  = lobster.biology.excretion_inorganic_fraction
    μ   = lobster.biology.zooplankton_excretion_rate
    αP  = lobster.biology.ammonia_fraction_of_exudate
    γ   = lobster.biology.phytoplankton_exudation_fraction
    
    @inbounds begin
        Z = fields.Z[i, j, k]
    end
    
    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    return (1 - αP) * γ * μP + (1 - αZ) * μ * Z
end

@inline biology_organic_carbon_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    lobster.biology.redfield_ratio * biology_organic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)

@inline function solid_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    αᴾ = lobster.biology.phytoplankton_solid_waste_fraction
    αᶻ = lobster.biology.zooplankton_assimilation_fraction
    mᴾ = lobster.biology.phytoplankton_mortality_rate
    mᶻ = lobster.biology.zooplankton_mortality_rate

    @inbounds begin
        P = fields.P[i, j, k]
        Z = fields.Z[i, j, k]
    end

    G = total_grazing(lobster, i, j, k, fields, auxiliary_fields)

    νP = mortality(lobster.biology.phytoplankton_mortality_formulation, P, mᴾ)

    return (1 - αᶻ) * G + αᴾ * νP + mᶻ * Z^2
end

@inline solid_carbon_waste(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields) =
    solid_waste(lobster, i, j, k, fields, auxiliary_fields) * lobster.biology.redfield_ratio

#### grazing
@inline function grazing(lobster::PHYTO_ZOO_BND, i, j, k, ::Val{:P}, fields, auxiliary_fields)
    kG = lobster.biology.grazing_half_saturation
    g  = lobster.biology.maximum_grazing_rate

    @inbounds begin
        Z = fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(lobster.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    L = concentration_limit(lobster.biology.grazing_concentration_formulation, food_concentration, kG)

    return g * p * L * P / (food_concentration + eps(food_concentration)) * Z
end

@inline function grazing(lobster::PHYTO_ZOO_BND, i, j, k, ::Union{Val{:sPOM}, Val{:sPON}, Val{:D}}, fields, auxiliary_fields)
    kG = lobster.biology.grazing_half_saturation
    g  = lobster.biology.maximum_grazing_rate

    @inbounds begin
        Z = fields.Z[i, j, k]
        P = fields.P[i, j, k]
        sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    end

    p = weighted_phytoplankton_preference(lobster.biology, P, sPOM)

    food_concentration = p * P + (1-p) * sPOM

    L = concentration_limit(lobster.biology.grazing_concentration_formulation, food_concentration, kG)

    return g * (1-p) * L * sPOM / (food_concentration + eps(food_concentration)) * Z
end

@inline grazing(lobster::PHYTO_ZOO_BND, i, j, k, ::Val{:sPOC}, fields, auxiliary_fields) =
    grazing(lobster, i, j, k, Val(:sPOM), fields, auxiliary_fields) * lobster.biology.redfield_ratio

@inline function weighted_phytoplankton_preference(biology::PhytoZoo, P, sPOM)
    p̃ = biology.preference_for_phytoplankton

    return p̃ * P / (p̃ * P + (1 - p̃) * sPOM + eps(0.0))
end

@inline function calcite_production(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    α  = lobster.biology.zooplankton_assimilation_fraction
    mᴾ = lobster.biology.phytoplankton_mortality_rate
    R  = lobster.biology.redfield_ratio 
    ρ  = lobster.biology.carbon_calcate_ratio
    η  = lobster.biology.zooplankton_gut_calcite_dissolution

    P = @inbounds fields.P[i, j, k]

    G = grazing(lobster, i, j, k, Val(:P), fields, auxiliary_fields)
    ν = mortality(lobster.biology.phytoplankton_mortality_formulation, P, mᴾ)

    return (G * (1 - η) + ν) * ρ * R
end

@inline function calcite_dissolution(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    R  = lobster.biology.redfield_ratio 
    ρ  = lobster.biology.carbon_calcate_ratio
    η  = lobster.biology.zooplankton_gut_calcite_dissolution

    G = grazing(lobster, i, j, k, Val(:P), fields, auxiliary_fields)

    return G * η * ρ * R
end

# when we have variable redfield particles the calcite production goes to bPOC, but this
# isn't implicitly captured in the bPOM when we have a fixed redfield, so we have to 
# assume instant dissolution and put it back in DIC (or we could choose to lose it),
# maybe that would be better?
# but when we don't this can't be implicitly captured so we put it all back in the DIC compartement
@inline function calcite_dissolution(lobster::BiologyNutrientDetritus{<:Any, <:PhytoZoo, <:Union{TwoParticleAndDissolved, Detritus, Nothing}}, i, j, k, fields, auxiliary_fields)
    R  = lobster.biology.redfield_ratio 
    ρ  = lobster.biology.carbon_calcate_ratio
    mᴾ = lobster.biology.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    G = grazing(lobster, i, j, k, Val(:P), fields, auxiliary_fields)
    ν = mortality(lobster.biology.phytoplankton_mortality_formulation, P, mᴾ)


    return (G + ν) * ρ * R
end

@inline function calcite_uptake(lobster::PHYTO_ZOO_BND, i, j, k, fields, auxiliary_fields)
    R  = lobster.biology.redfield_ratio 
    ρ  = lobster.biology.carbon_calcate_ratio

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    return 2 * ρ * μP * R
end