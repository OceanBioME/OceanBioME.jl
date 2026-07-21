using Oceananigans.Fields: ZeroField

"""
    PhytoZoo([FT = Float64;] kwargs...)
    PhytoZoo(grid; phytoplankton_sinking_speed = 0, zooplankton_sinking_speed = 0, open_bottom = true, kwargs...)

A phytoplankton-and-zooplankton plankton component for the `plankton` slot of a
[`NutrientsPlanktonDetritus`](@ref) model, as used by the [`LOBSTER`](@ref) and
[`NPZD`](@ref) presets. It adds phytoplankton (`P`) and zooplankton (`Z`) tracers (plus a
temperature tracer `T` when a `temperature_coefficient` is supplied). Phytoplankton grow with light-
and nutrient-limited uptake, are grazed by zooplankton, and both pools produce dissolved/solid waste.

Use the `grid` method when the plankton should sink; it configures the sinking-velocity fields from
`phytoplankton_sinking_speed`/`zooplankton_sinking_speed` and forwards the rest to the base
constructor.

Keyword Arguments
=================

- `nutrient_half_saturations`: a `NamedTuple` of half-saturation constants keyed by limiting nutrient
  (e.g. `nitrate`, `ammonia`, `iron`); its keys set which nutrients limit growth
- `light_limitation`: the light-limitation functional form (`MondoLightLimitation()` or
  `AnalyticalLightLimitation()`), with `light_half_saturation` (W/m²)
- `phytoplankton_maximum_growth_rate`, `phytoplankton_exudation_fraction`,
  `phytoplankton_mortality_rate`, `phytoplankton_solid_waste_fraction`: phytoplankton growth and loss
  parameters
- `maximum_grazing_rate`, `grazing_half_saturation`, `preference_for_phytoplankton`,
  `zooplankton_assimilation_fraction`, `zooplankton_mortality_rate`, `zooplankton_excretion_rate`:
  zooplankton grazing and loss parameters
- `temperature_coefficient`: optional Q₁₀ temperature dependence (off by default; adds a `T` tracer)
- `carbon_ratio`, `iron_ratio`, `phosphate_ratio`, `rain_ratio`, `chlorophyll_ratio`: elemental and
  chlorophyll ratios relative to nitrogen
- see the constructor definition for the full list of tunable parameters and their default values/units
"""
struct PhytoZoo{LN, EN, TT, FT, LL, TC, PS, ZS} <: AbstractPlankton{LN}
          nutrient_half_saturations :: NamedTuple{LN, TT}

         nitrate_ammonia_inhibition :: FT
              light_half_saturation :: FT
                   light_limitation :: LL

  phytoplankton_maximum_growth_rate :: FT
   phytoplankton_exudation_fraction :: FT 
        ammonia_fraction_of_exudate :: FT 

            temperature_coefficient :: TC

       phytoplankton_mortality_rate :: FT
         zooplankton_mortality_rate :: FT
         zooplankton_excretion_rate :: FT

 phytoplankton_solid_waste_fraction :: FT
       excretion_inorganic_fraction :: FT 
       preference_for_phytoplankton :: FT 
               maximum_grazing_rate :: FT
            grazing_half_saturation :: FT
  zooplankton_assimilation_fraction :: FT
        edible_fraction_of_detritus :: FT
               edible_detritus_name :: Val{EN}

# In the original this fraction of calcite is dissolved in zooplankton, but without explicitly tracking calcite we can't close the
# budget because the calcite to carbon ratio in the detritus will vary depending on where the waste comes from
#    zooplankton_calcite_dissolution :: FT = 0.3

                       carbon_ratio :: FT
                         iron_ratio :: FT
                    phosphate_ratio :: FT

                         rain_ratio :: FT

                  chlorophyll_ratio :: FT 

     phytoplankton_sinking_velocity :: PS
       zooplankton_sinking_velocity :: ZS
end

function PhytoZoo(FT = Float64;
                  nutrient_half_saturations = 
                      (nitrate = 0.7,                              # mmol N/m³
                       ammonia = 0.001,                            # mmol N/m³
                       iron = 2e-4),                               # mmol Fe/m³ - estimated from fitting OSP

                  nitrate_ammonia_inhibition = 3.0,       
                  light_half_saturation = 33.0,                  # W/m² 
                  light_limitation = MondoLightLimitation(),

                  phytoplankton_maximum_growth_rate = 2.42e-5,    # 1/s (double documented due to change to nitrogen limitation)

                  phytoplankton_exudation_fraction = 0.05,
                  ammonia_fraction_of_exudate = 0.75,

                  temperature_coefficient = nothing,          # Q10 factor, off by default 1.88 option from Kuhn, 2015

                  phytoplankton_mortality_rate = 5.8e-7,          # 1/s/mmol N/m³
                  zooplankton_mortality_rate = 2.31e-6,           # 1/s/mmol N/m³
                  zooplankton_excretion_rate = 5.8e-7,            # 1/s

                  phytoplankton_solid_waste_fraction = 1.0,
                  excretion_inorganic_fraction = 0.5,
                  preference_for_phytoplankton = 0.5,
                  maximum_grazing_rate = 9.26e-6,                 # 1/s
                  grazing_half_saturation = 1.0,                  # mmol N/m³
                  zooplankton_assimilation_fraction = 0.7,
                  edible_fraction_of_detritus = 0.5,
                  edible_detritus_name = :sPOM,                   # if multiple detritus classes are present
                                                                  # this is edible in full, if only one then the fraction is used
                  carbon_ratio = 6.56,                            # mol C / mol N
                  iron_ratio = 4.6375e-5,                         # mol Fe / mol N - from PISCES optimal ratio
                  phosphate_ratio = 1/16,                         # mol P / mol N

                  rain_ratio = 0.1,                               # mol CaCO₃ / mol C

                  chlorophyll_ratio = 1.31,                       # g Chl/mol N

                  phytoplankton_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()),
                  zooplankton_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()))

     limiting_nutrients = keys(nutrient_half_saturations)

     nutrient_half_saturations = NamedTuple{limiting_nutrients}(map(val->convert(FT, val), values(nutrient_half_saturations)))

     temperature_coefficient = isnothing(temperature_coefficient) ? nothing : convert(FT, temperature_coefficient)

     return PhytoZoo(nutrient_half_saturations,
                     convert(FT, nitrate_ammonia_inhibition),
                     convert(FT, light_half_saturation),
                     light_limitation,
                     convert(FT, phytoplankton_maximum_growth_rate),
                     convert(FT, phytoplankton_exudation_fraction),
                     convert(FT, ammonia_fraction_of_exudate),
                     temperature_coefficient,
                     convert(FT, phytoplankton_mortality_rate),
                     convert(FT, zooplankton_mortality_rate),
                     convert(FT, zooplankton_excretion_rate),
                     convert(FT, phytoplankton_solid_waste_fraction),
                     convert(FT, excretion_inorganic_fraction),
                     convert(FT, preference_for_phytoplankton),
                     convert(FT, maximum_grazing_rate),
                     convert(FT, grazing_half_saturation),
                     convert(FT, zooplankton_assimilation_fraction),
                     convert(FT, edible_fraction_of_detritus),
                     Val(edible_detritus_name),
                     convert(FT, carbon_ratio),
                     convert(FT, iron_ratio),
                     convert(FT, phosphate_ratio),
                     convert(FT, rain_ratio),
                     convert(FT, chlorophyll_ratio),
                     phytoplankton_sinking_velocity,
                     zooplankton_sinking_velocity)
end

function PhytoZoo(grid::AbstractGrid{FT};
                  phytoplankton_sinking_speed = zero(FT),
                  zooplankton_sinking_speed = zero(FT),
                  open_bottom = true,
                  kwargs...) where FT

    sinking_velocities = setup_velocity_fields((P = phytoplankton_sinking_speed, 
                                                Z = zooplankton_sinking_speed), 
                                                grid, open_bottom; three_D = true)

    phytoplankton_sinking_velocity = sinking_velocities.P
    zooplankton_sinking_velocity = sinking_velocities.Z

    return PhytoZoo(FT; phytoplankton_sinking_velocity,
                        zooplankton_sinking_velocity,
                        kwargs...)
end

const PhytoZoo_NPD{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:PhytoZoo}

required_biogeochemical_tracers(pz::PhytoZoo) = ((:P, :Z)..., (isnothing(pz.temperature_coefficient) ? () : (:T, ))...)
required_biogeochemical_auxiliary_fields(::PhytoZoo) = (:PAR, )

biogeochemical_drift_velocity(bgc::PhytoZoo_NPD, ::Val{:P}) = 
    bgc.plankton.phytoplankton_sinking_velocity

biogeochemical_drift_velocity(bgc::PhytoZoo_NPD, ::Val{:Z}) = 
    bgc.plankton.zooplankton_sinking_velocity

@inline nitrogen_ratio(::PhytoZoo, ::NPD{FT}) where FT = one(FT)
@inline phosphate_ratio(plankton::PhytoZoo, ::NPD{FT}) where FT = plankton.phosphate_ratio
@inline carbon_ratio(plankton::PhytoZoo, ::NPD{FT}) where FT = plankton.carbon_ratio
@inline iron_ratio(plankton::PhytoZoo, ::NPD{FT}) where FT = plankton.iron_ratio
@inline calcite_rain_ratio(plankton::PhytoZoo, ::NPD{FT}) where FT = plankton.rain_ratio
@inline chlorophyll(plankton::PhytoZoo, model) = plankton.chlorophyll_ratio * model.tracers.P

@inline limiting_nutrients(::PhytoZoo{LN}) where LN = LN
@inline edible_detritus_name(::PhytoZoo{<:Any, EN}) where EN = EN

@inline nutrient_half_saturations(plankton::PhytoZoo, ::Val{:N})   = plankton.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(plankton::PhytoZoo, ::Val{:NO₃}) = plankton.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(plankton::PhytoZoo, ::Val{:NH₄}) = plankton.nutrient_half_saturations.ammonia
@inline nutrient_half_saturations(plankton::PhytoZoo, ::Val{:PO₄}) = plankton.nutrient_half_saturations.phosphate
@inline nutrient_half_saturations(plankton::PhytoZoo, ::Val{:Fe})  = plankton.nutrient_half_saturations.iron

@inline function (bgc::PhytoZoo_NPD)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    γ = bgc.plankton.phytoplankton_exudation_fraction
    m = bgc.plankton.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    μP = phytoplankton_growth(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

    Gp = grazing(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)

    return (1 - γ) * μP - Gp - m * P^2
end 

@inline function (bgc::PhytoZoo_NPD)(i, j, k, grid, ::Val{:Z}, clock, fields, auxiliary_fields)
    α = bgc.plankton.zooplankton_assimilation_fraction
    m = bgc.plankton.zooplankton_mortality_rate
    μ = bgc.plankton.zooplankton_excretion_rate

    Z = @inbounds fields.Z[i, j, k]
    G = total_grazing(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

    return α * G - m * Z^2 - μ * Z
end 

# growth and nutrient uptake
@inline function phytoplankton_growth(i, j, k, grid, plankton::PhytoZoo, bgc, fields, auxiliary_fields)
    kPAR = plankton.light_half_saturation
    μ₀   = plankton.phytoplankton_maximum_growth_rate
    Q10  = plankton.temperature_coefficient

    @inbounds begin
        PAR = auxiliary_fields.PAR[i, j, k]
        P   = fields.P[i, j, k]
    end

    Lₙ = nutrient_limitation(i, j, k, grid, bgc.nutrients, plankton, bgc, fields, auxiliary_fields)
    Lₗ = light_limitation(plankton.light_limitation, PAR, kPAR)
    Lₜ = temperature_limitation(Q10, i, j, k, fields)

    return μ₀ * Lₗ * Lₙ * Lₜ * P
end

@inline temperature_limitation(::Nothing, i, j, k, fields) = one(eltype(fields.P))

@inline function temperature_limitation(Q10::FT, i, j, k, fields) where FT
    T = @inbounds fields.T[i, j, k]

    return Q10 ^ (T / convert(FT, 10))
end

@inline function nitrogen_limitation(i, j, k, grid, 
                                     ::NitrateAmmonia, 
                                     plankton::PhytoZoo,
                                     ::NutrientsPlanktonDetritus{FT},
                                     fields, auxiliary_fields) where FT

    kNO₃ = nutrient_half_saturations(plankton, Val(:NO₃)) 
    kNH₄ = nutrient_half_saturations(plankton, Val(:NH₄))
    ψ    = plankton.nitrate_ammonia_inhibition

    @inbounds begin 
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    return (nitrate_limitation + ammonia_limitation) / 2
end

@inline nutrient_uptake(i, j, k, grid,
                        plankton::PhytoZoo, bgc,
                        fields, auxiliary_fields) = 
    phytoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

@inline function nutrient_uptake(i, j, k, grid, ::Val{:NO₃}, plankton::PhytoZoo, bgc::NPD{FT}, fields, auxiliary_fields) where FT
    μ = phytoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    kNO₃ = nutrient_half_saturations(plankton, Val(:NO₃)) 
    kNH₄ = nutrient_half_saturations(plankton, Val(:NH₄))
    ψ    = plankton.nitrate_ammonia_inhibition

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    return μ * nitrate_limitation / (nitrate_limitation + ammonia_limitation + eps(zero(FT)))
end

@inline function nutrient_uptake(i, j, k, grid, ::Val{:NH₄}, plankton::PhytoZoo, bgc::NPD{FT}, fields, auxiliary_fields) where FT
    μ = phytoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    kNO₃ = nutrient_half_saturations(plankton, Val(:NO₃)) 
    kNH₄ = nutrient_half_saturations(plankton, Val(:NH₄))
    ψ    = plankton.nitrate_ammonia_inhibition
    
    @inbounds begin
        NO₃ = fields.NO₃[i, j, k]
        NH₄ = fields.NH₄[i, j, k]
    end

    nitrate_limitation = NO₃ * exp(-ψ * NH₄) / (NO₃ + kNO₃)
    ammonia_limitation = max(0, NH₄ / (kNH₄ + NH₄))

    return μ * ammonia_limitation / (nitrate_limitation + ammonia_limitation + eps(zero(FT)))
end

# grazing
@inline function total_grazing(i, j, k, grid, plankton::PhytoZoo, bgc, fields, auxiliary_fields)
    kG = plankton.grazing_half_saturation
    g  = plankton.maximum_grazing_rate

    @inbounds begin
        Z    = fields.Z[i, j, k]
        P    = fields.P[i, j, k]
        ePOM = edible_particulate_organic_matter(i, j, k, grid, bgc.detritus, bgc.plankton, bgc, fields)
    end

    p = weighted_phytoplankton_preference(bgc.plankton, P, ePOM)

    food_concentration = p * P + (1-p) * ePOM

    L = food_concentration^2 / (food_concentration^2 + kG^2)

    return g * L * Z
end

@inline function grazing(i, j, k, grid, ::Val{:P}, plankton::PhytoZoo, bgc::NPD, fields, auxiliary_fields)
    kG = plankton.grazing_half_saturation
    g  = plankton.maximum_grazing_rate

    @inbounds begin
        Z    = fields.Z[i, j, k]
        P    = fields.P[i, j, k]
        ePOM = edible_particulate_organic_matter(i, j, k, grid, bgc.detritus, bgc.plankton, bgc, fields)
    end

    p = weighted_phytoplankton_preference(bgc.plankton, P, ePOM)

    food_concentration = p * P + (1-p) * ePOM

    L = food_concentration^2 / (food_concentration^2 + kG^2)

    return g * p * L * P / (food_concentration + eps(zero(P))) * Z
end

@inline function grazing(i, j, k, grid, ::Union{Val{:sPOM}, Val{:sPON}, Val{:sPOP}, Val{:D}, Val{:POM}, Val{:POP}}, plankton::PhytoZoo, bgc::NPD, fields, auxiliary_fields)
    kG = plankton.grazing_half_saturation
    g  = plankton.maximum_grazing_rate

    @inbounds begin
        Z    = fields.Z[i, j, k]
        P    = fields.P[i, j, k]
        ePOM = edible_particulate_organic_matter(i, j, k, grid, bgc.detritus, bgc.plankton, bgc, fields)
    end

    p = weighted_phytoplankton_preference(bgc.plankton, P, ePOM)

    food_concentration = p * P + (1-p) * ePOM

    L = food_concentration^2 / (food_concentration^2 + kG^2)

    return g * (1-p) * L * ePOM / (food_concentration + eps(zero(P))) * Z
end

@inline grazing(i, j, k, grid, ::Val{:sPOC}, plankton::PhytoZoo, bgc::NPD, fields, auxiliary_fields) =
    carbon_ratio(i, j, k, grid, plankton, bgc, fields) * 
    grazing(i, j, k, grid, Val(:sPON), plankton, bgc, fields, auxiliary_fields)

@inline function weighted_phytoplankton_preference(plankton::PhytoZoo, P::FT, sPOM) where FT
    p̃ = plankton.preference_for_phytoplankton

    return p̃ * P / (p̃ * P + (1 - p̃) * sPOM + eps(zero(FT)))
end

# assumptions about particle edibility
@inline edible_particulate_organic_matter(i, j, k, grid, ::InstantRemineralisationDetritus, plankton::PhytoZoo, bgc::NPD{FT}, fields) where FT = 
    zero(FT)
    
@inline edible_particulate_organic_matter(i, j, k, grid, ::Detritus, plankton::PhytoZoo, bgc, fields) = 
    @inbounds plankton.edible_fraction_of_detritus * fields.D[i, j, k]

@inline edible_particulate_organic_matter(i, j, k, grid, ::DissolvedParticulate{<:Any, 1, <:Any, PN}, plankton::PhytoZoo, bgc, fields) where PN = 
    @inbounds plankton.edible_fraction_of_detritus * getproperty(fields, PN[1])[i, j, k]

@inline edible_particulate_organic_matter(i, j, k, grid, ::DissolvedParticulate, plankton::PhytoZoo, bgc, fields) =
    @inbounds getproperty(fields, edible_detritus_name(plankton))[i, j, k] 

@inline edible_particulate_organic_matter(i, j, k, grid, ::CarbonNitrogenDissolvedParticulate, plankton::PhytoZoo, bgc, fields) = 
    @inbounds fields.sPON[i, j, k] 

# waste routing
@inline function solid_waste(i, j, k, grid, plankton::PhytoZoo, bgc, fields, auxiliary_fields)
    αᴾ = plankton.phytoplankton_solid_waste_fraction
    αᶻ = plankton.zooplankton_assimilation_fraction
    mᴾ = plankton.phytoplankton_mortality_rate
    mᶻ = plankton.zooplankton_mortality_rate

    @inbounds begin
        P = fields.P[i, j, k]
        Z = fields.Z[i, j, k]
    end

    G = total_grazing(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    return (1 - αᶻ) * G + αᴾ * mᴾ * P^2 + mᶻ * Z^2
end

@inline function dissolved_waste(i, j, k, grid, plankton::PhytoZoo, bgc, fields, auxiliary_fields)
    αᶻ  = plankton.excretion_inorganic_fraction
    μ   = plankton.zooplankton_excretion_rate
    αᴾ  = plankton.ammonia_fraction_of_exudate
    γ   = plankton.phytoplankton_exudation_fraction
    
    @inbounds begin
        Z = fields.Z[i, j, k]
    end
    
    μP = phytoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    return (1 - αᴾ) * γ * μP + (1 - αᶻ) * μ * Z
end

@inline function inorganic_waste(i, j, k, grid, plankton::PhytoZoo, bgc, fields, auxiliary_fields)
    αᴾ = plankton.phytoplankton_solid_waste_fraction
    αᶻ = plankton.excretion_inorganic_fraction
    μ  = plankton.zooplankton_excretion_rate
    mᴾ = plankton.phytoplankton_mortality_rate
    αᵃ = plankton.ammonia_fraction_of_exudate
    γ  = plankton.phytoplankton_exudation_fraction

    @inbounds begin
        P = fields.P[i, j, k]
        Z = fields.Z[i, j, k]
    end

    μP = phytoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    return αᶻ * μ * Z + (1 - αᴾ) * mᴾ * P^2 + αᵃ * γ * μP
end

# admin
Adapt.adapt_structure(to, pz::PhytoZoo) = 
    PhytoZoo(adapt(to, pz.nutrient_half_saturations),
             adapt(to, pz.nitrate_ammonia_inhibition),
             adapt(to, pz.light_half_saturation),
             adapt(to, pz.light_limitation),
             adapt(to, pz.phytoplankton_maximum_growth_rate),
             adapt(to, pz.phytoplankton_exudation_fraction),
             adapt(to, pz.ammonia_fraction_of_exudate),
             adapt(to, pz.temperature_coefficient),
             adapt(to, pz.phytoplankton_mortality_rate),
             adapt(to, pz.zooplankton_mortality_rate),
             adapt(to, pz.zooplankton_excretion_rate),
             adapt(to, pz.phytoplankton_solid_waste_fraction),
             adapt(to, pz.excretion_inorganic_fraction),
             adapt(to, pz.preference_for_phytoplankton),
             adapt(to, pz.maximum_grazing_rate),
             adapt(to, pz.grazing_half_saturation),
             adapt(to, pz.zooplankton_assimilation_fraction),
             adapt(to, pz.edible_fraction_of_detritus),
             adapt(to, pz.edible_detritus_name),
             adapt(to, pz.carbon_ratio),
             adapt(to, pz.iron_ratio),
             adapt(to, pz.phosphate_ratio),
             adapt(to, pz.rain_ratio),
             adapt(to, pz.chlorophyll_ratio),
             adapt(to, pz.phytoplankton_sinking_velocity),
             adapt(to, pz.zooplankton_sinking_velocity))

Base.summary(::PhytoZoo) = "PhytoZoo (:P, :Z)"

function Base.show(io::IO, pz::PhytoZoo)
    msg = summary(pz) * "\n"

    msg *= "├── C:N:P:Fe = $(pz.carbon_ratio):1:$(pz.phosphate_ratio):$(pz.iron_ratio)\n"
    msg *= "└── Limiting nutrients: $(limiting_nutrients(pz))"

    print(io, msg)

    return nothing
end
