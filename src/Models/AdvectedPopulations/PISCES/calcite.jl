"""
    Calcite

Stores the parameter values for calcite (`CaCO₃`) evolution.

Keyword Arguments
=================
- `base_rain_ratio`: the base fraction of Coccolithophores
- `base_dissolution_rate`: base rate of calcite dissolution (1/s)
- `dissolution_exponent`: exponent of calcite excess for dissolution rate 

"""
@kwdef struct Calcite{FT}
          base_rain_ratio :: FT = 0.3         # 
    base_dissolution_rate :: FT = 0.197 / day # 1 / s
     dissolution_exponent :: FT = 1.0         # 
end

@inline function (calcite::Calcite)(::Val{:CaCO₃}, bgc,
                                    x, y, z, t,
                                    P, D, Z, M, 
                                    PChl, DChl, PFe, DFe, DSi,
                                    DOC, POC, GOC, 
                                    SFe, BFe, PSi, 
                                    NO₃, NH₄, PO₄, Fe, Si, 
                                    CaCO₃, DIC, Alk, 
                                    O₂, T, S,
                                    zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

    production = calcite_production(calcite, bgc, z, P, D, PChl, PFe, Z, M, POC, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)
    
    dissolution = calcite_dissolution(calcite, CaCO₃, Ω)
    
    return production - dissolution
end

@inline function calcite_production(calcite, bgc, z, P, D, PChl, PFe, Z, M, POC, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)
    R = rain_ratio(calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    microzooplankton = specific_calcite_grazing_loss(bgc.microzooplankton, P, D, Z, POC, T) * Z
    mesozooplankton  = specific_calcite_grazing_loss(bgc.mesozooplankton, P, D, Z, POC, T) * M

    linear_mortality, quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)

    total_mortality = 0.5 * (linear_mortality + quadratic_mortality)

    return R * (microzooplankton + mesozooplankton + total_mortality)
end

# should this be in the particles thing?
@inline function rain_ratio(calcite::Calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)
    r = calcite.base_rain_ratio

    # assuming this is a type in Aumont 2015 based on Aumont 2005
    _, _, LPO₄, LN = bgc.nanophytoplankton.nutrient_limitation(bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    LFe = Fe / (Fe + 0.05)

    # from NEMO source code +/- presumably a typo replacing kPO4 with kNH4 which aren't even same units
    L_CaCO₃ = min(LN, LPO₄, LFe)
    
    phytoplankton_concentration_factor = max(1, P / 2)

    low_light_factor = max(0, PAR - 1) / (4 + PAR)
    high_light_factor = 30 / (30 + PAR)

    low_temperature_factor = max(0, T / (T + 0.1)) # modified from origional as it goes negative and does not achieve goal otherwise
    high_temperature_factor = 1 + exp(-(T - 10)^2 / 25)

    depth_factor = min(1, -50/zₘₓₗ)

    return r * L_CaCO₃ * phytoplankton_concentration_factor * low_light_factor * high_light_factor * low_temperature_factor * high_temperature_factor * depth_factor
end

@inline function calcite_dissolution(calcite, CaCO₃, Ω)
    λ   = calcite.base_dissolution_rate
    nca = calcite.dissolution_exponent

    ΔCaCO₃ = max(0, 1 - Ω)

    return λ * ΔCaCO₃ ^ nca * CaCO₃
end
