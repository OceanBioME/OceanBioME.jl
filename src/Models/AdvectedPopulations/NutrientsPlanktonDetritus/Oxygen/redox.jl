struct RedoxOxygen{FT, SF}
              nitrate_carbon_to_oxygen :: FT  # mol C/mol O₂, nitrate-based new production; MARBL parm_Red_D_C_O2
     remineralisation_carbon_to_oxygen :: FT  # mol C/mol O₂, remin & ammonium-based production; MARBL parm_Remin_D_C_O2
           diazotroph_carbon_to_oxygen :: FT  # mol C/mol O₂, N-fixation-based production; MARBL parm_Red_D_C_O2_diaz
    denitrification_carbon_to_nitrogen :: FT  # mol C/mol N, organic C remineralised on NO₃; MARBL denitrif_C_N

                    nitrification_rate :: FT  # 1/s, NH₄→NO₃ in the dark; MARBL parm_kappa_nitrif
               nitrification_par_limit :: FT  # W/m², nitrification stops above this; MARBL parm_nitrif_par_lim
                        oxygen_minimum :: FT  # mmol/m³; MARBL parm_o2_min
                  oxygen_minimum_range :: FT  # mmol/m³, width of the taper; MARBL parm_o2_min_delta

           fastest_depletion_timescale :: FT
       oxygen_consumption_scale_factor :: SF  # prescribed forcing multiplying O₂ consumption; MARBL o2_consumption_scalef
end

"""
    RedoxOxygen([FT = Float64;] kwargs...)

An oxygen component for the `oxygen` slot of a [`NutrientsPlanktonDetritus`](@ref) model which also
carries nitrogen redox chemistry. It adds the oxygen tracer (`O₂`) and, unlike [`Oxygen`](@ref),
feeds back on nitrogen: ammonium is nitrified to nitrate in dim water, and in suboxic water nitrate is
denitrified to N₂ as organic carbon is remineralised.

Oxygen production is split by nitrogen source — nitrate-based new production releases more oxygen per
carbon than ammonium-based regenerated production, and nitrogen fixation more still — so it requires a
plankton component that reports uptake per source.

Keyword Arguments
=================

- `nitrate_carbon_to_oxygen`, `remineralisation_carbon_to_oxygen`, `diazotroph_carbon_to_oxygen`:
  carbon-to-oxygen ratios for each production and remineralisation pathway (mol C/mol O₂)
- `denitrification_carbon_to_nitrogen`: carbon remineralised per nitrogen denitrified (mol C/mol N)
- `nitrification_rate_per_day`, `nitrification_par_limit`: the dark nitrification rate (1/day) and the
  light level above which it stops (W/m²)
- `oxygen_minimum`, `oxygen_minimum_range`: the oxygen concentration below which respiration gives way
  to denitrification, and the width of that transition (mmol/m³)
- `oxygen_consumption_scale_factor`: a prescribed scalar or 3-D `Field` multiplying oxygen consumption
  (default `1`, i.e. no scaling)
"""
function RedoxOxygen(FT = Float64;
               nitrate_carbon_to_oxygen           = 117/170,
               remineralisation_carbon_to_oxygen  = 117/138,
               diazotroph_carbon_to_oxygen        = 117/150,
               denitrification_carbon_to_nitrogen = 117/136,
               nitrification_rate_per_day         = 0.06,   # 1/day
               nitrification_par_limit            = 1.0,    # W/m²
               oxygen_minimum                     = 5.0,    # mmol/m³
               oxygen_minimum_range               = 5.0,    # mmol/m³
               fastest_depletion_timescale        = 10days,
               oxygen_consumption_scale_factor    = 1.0)

    scale_factor = oxygen_consumption_scale_factor isa Number ?
                   convert(FT, oxygen_consumption_scale_factor) : oxygen_consumption_scale_factor

    return RedoxOxygen{FT, typeof(scale_factor)}(convert(FT, nitrate_carbon_to_oxygen),
                                           convert(FT, remineralisation_carbon_to_oxygen),
                                           convert(FT, diazotroph_carbon_to_oxygen),
                                           convert(FT, denitrification_carbon_to_nitrogen),
                                           convert(FT, nitrification_rate_per_day / day),
                                           convert(FT, nitrification_par_limit),
                                           convert(FT, oxygen_minimum),
                                           convert(FT, oxygen_minimum_range),
                                           convert(FT, fastest_depletion_timescale),
                                           scale_factor)
end

required_biogeochemical_tracers(::RedoxOxygen) = (:O₂, )
required_biogeochemical_auxiliary_fields(::RedoxOxygen) = (:PAR, )

# needed for sediment scaling
biogeochemical_auxiliary_fields(ox::RedoxOxygen) =
    (O₂_consumption_scale = ox.oxygen_consumption_scale_factor isa AbstractField ?
                                ox.oxygen_consumption_scale_factor : 
                                ConstantField(o.oxygen_consumption_scale_factor), )

const RedoxOxygen_NPD = NutrientsPlanktonDetritus{<:Any, <:Nutrients{<:NitrateAmmonia}, <:Any,
                                                  <:AbstractRefractoryDissolvedDetritus, <:Any, 
                                                  <:RedoxOxygen}

#####
##### Pure rate helpers
#####

@subcolumn_average @inline function nitrification_taper(PAR_top, PAR_bot, PAR_nit)
    PAR_bot = max(PAR_bot, eps(zero(PAR_bot)))

    partial = log(PAR_bot / PAR_nit) / log(PAR_bot / PAR_top)

    return ifelse(PAR_bot >= I_nit, zero(partial),
                  ifelse(PAR_top <= PAR_nit, one(partial), partial))
end

@inline suboxic_fraction(O₂, o2_min, Δ) = clamp(((o2_min + Δ) - O₂) / Δ, zero(O₂), one(O₂))
@inline aerobic_fraction(O₂, o2_min, Δ) = clamp((O₂ - o2_min) / Δ, zero(O₂), one(O₂))

#####
##### Coupling functions
#####
@inline function water_column_nitrification(i, j, k, grid, bgc, fields, aux)
    ox = bgc.oxygen

    nh = @inbounds fields.NH₄[i, j, k]

    NH₄ = max(nh, zero(nh))

    I_bot = interface_par(i, j, k,     grid, aux)    # the deeper face
    I_top = interface_par(i, j, k + 1, grid, aux)    # the shallower face
    I_nit = ox.nitrification_par_limit

    taper = nitrification_taper(I_top, I_bot, I_nit) # subcolumn average if available

    return ox.nitrification_rate * NH₄ * taper
end

@inline function water_column_denitrification(i, j, k, grid, bgc, fields, aux)
    ox = bgc.oxygen

    @inbounds begin
        O₂  = fields.O₂[i, j, k]
        NO₃ = fields.NO₃[i, j, k]
    end

    ω     = suboxic_fraction(O₂, ox.oxygen_minimum, ox.oxygen_minimum_range)
    remin = inorganic_carbon_waste(i, j, k, grid, bgc.detritus, bgc, fields, aux)

    raw = ω * remin / ox.denitrification_carbon_to_nitrogen

    return min(raw, max(NO₃, zero(NO₃)) / ox.fastest_depletion_timescale)
end

# this gets scaled 
@inline function oxygen_consumption(i, j, k, grid, bgc, fields, aux)
    o = bgc.oxygen

    @inbounds O₂ = fields.O₂[i, j, k]

    η = aerobic_fraction(O₂, o.oxygen_minimum, o.oxygen_minimum_range)

    remin = inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) +
            inorganic_carbon_waste(i, j, k, grid, bgc.detritus, bgc, fields, aux)

    nitrif = water_column_nitrification(i, j, k, grid, bgc, fields, aux)

    return η * (remin / o.remineralisation_carbon_to_oxygen + 2 * nitrif)
end

@inline get_oxygen_consumption_scale(i, j, k, grid, Φ) = @inbounds Φ[i, j, k]
@inline get_oxygen_consumption_scale(i, j, k, grid, Φ::Number) = Φ

@inline (bgc::RedoxOxygen_NPD)(i, j, k, grid, ::Val{:O₂}, clock, fields, aux) =
    oxygen_production_total(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    get_oxygen_consumption_scale(i, j, k, grid, bgc.oxygen.oxygen_consumption_scale_factor) *
    oxygen_consumption(i, j, k, grid, bgc, fields, aux)

@inline (bgc::RedoxOxygen_NPD)(i, j, k, grid, ::Val{:NO₃}, clock, fields, aux) =
    water_column_nitrification(i, j, k, grid, bgc, fields, aux) -
    water_column_denitrification(i, j, k, grid, bgc, fields, aux) -
    nutrient_uptake(i, j, k, grid, Val(:NO₃), bgc.plankton, bgc, fields, aux)

@inline (bgc::RedoxOxygen_NPD)(i, j, k, grid, ::Val{:NH₄}, clock, fields, aux) =
    inorganic_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) +
    inorganic_nitrogen_waste(i, j, k, grid, bgc.detritus, bgc, fields, aux) -
    nutrient_uptake(i, j, k, grid, Val(:NH₄), bgc.plankton, bgc, fields, aux) -
    water_column_nitrification(i, j, k, grid, bgc, fields, aux)

@inline ballast_oxygen_scale(i, j, k, grid, oxygen::RedoxOxygen, ballast, fields) =
    @inbounds oxygen_scale_factor(fields.O₂[i, j, k], ballast)
