#####
##### Mineral-protection (ballast) sinking parameters
#####
##### Armstrong et al. (2000): a fixed mass of organic carbon rides down protected by each unit of
##### sinking mineral, and so travels with the mineral's much longer dissolution length scale rather than
##### its own. Each mineral is split into a "soft" class, which dissolves over a depth- and oxygen-
##### dependent length, and a "hard" class, which is essentially inert.
#####
##### Length scales are in metres.
#####

struct Ballast{FT, SL, D, SF}
    # soft-class base dissolution lengths (m)
    particulate_organic_dissolution_length :: FT
                calcite_dissolution_length :: FT
                   opal_dissolution_length :: FT
                   dust_dissolution_length :: FT

    # hard-class dissolution lengths (m); independent of depth and oxygen
         hard_dissolution_length :: FT
    dust_hard_dissolution_length :: FT

    # the fraction of production routed to the hard sub-class
    calcite_hard_fraction :: FT
       opal_hard_fraction :: FT
       dust_hard_fraction :: FT

    # hard organic carbon carried per unit mineral
    calcite_ballast_ratio :: FT
       opal_ballast_ratio :: FT
       dust_ballast_ratio :: FT

    # molar masses (g/mol)
    particulate_organic_carbon_mass :: FT
                       calcite_mass :: FT
                          opal_mass :: FT

    # the depth scale-length factor, piecewise linear through these points. Evaluated at the cell's
    # LOWER FACE rather than its centre — a discretisation choice the parameter values are tuned around.
    scale_length_depths :: SL   # m
    scale_length_values :: SL

    # low-oxygen stretch: slower remineralisation in an oxygen minimum zone
    oxygen_scaling_maximum :: FT  # mmol/m³, above which there is no stretch
    oxygen_scaling_minimum :: FT  # mmol/m³, below which the stretch is full
     oxygen_scaling_factor :: FT

    # particulate iron
        iron_desorption_rate :: FT  # per METRE — it multiplies a flux to give a volumetric rate
    iron_to_carbon_reference :: FT  # the fallback ratio where no organic flux comes in
                dust_to_iron :: FT  # mmol Fe per g of dust dissolved

    # boundary conditions and prescribed forcings
        surface_dust_flux :: D   # g/m²/s; scalar or field. Zero makes the whole dust pathway inert.
    sedimentary_iron_flux :: SF  # mmol Fe m⁻² s⁻¹; scalar or 3-D field
              open_bottom :: Bool  # true ⇒ the floor flux leaves the domain for a sediment to catch;
                                   # false ⇒ the bottom cell remineralises it, closing the column exactly
end

"""
    Ballast([FT = Float64;] kwargs...)

Mineral-protection sinking parameters for [`DissolvedRefractoryImplicitCNP`](@ref), after Armstrong et
al. (2000). All length scales are in metres.

Keyword Arguments
=================

- `*_dissolution_length`: the soft-class dissolution length of each mineral and of organic carbon (m)
- `hard_dissolution_length`, `dust_hard_dissolution_length`: the hard-class lengths (m)
- `*_hard_fraction`: the share of each mineral's production routed to the hard sub-class
- `*_hard_carbon_ratio`, `*_mass`: set how much organic carbon each unit of mineral protects
- `scale_length_depths`, `scale_length_values`: the piecewise-linear depth stretch of the soft lengths
- `oxygen_scaling_*`: the low-oxygen stretch
- `iron_desorption_rate`, `iron_to_carbon_reference`, `dust_to_iron`: the particulate iron budget
- `surface_dust_flux`, `sedimentary_iron_flux`: prescribed forcings (scalars or fields)
- `open_bottom`: whether the floor flux leaves the domain
"""
function Ballast(FT = Float64;
                 particulate_organic_dissolution_length = 100.0,   # m
                 calcite_dissolution_length             = 500.0,   # m
                 opal_dissolution_length                = 650.0,   # m
                 dust_dissolution_length                = 400.0,   # m
                 hard_dissolution_length                = 4.0e4,   # m
                 dust_hard_dissolution_length           = 1.2e6,   # m

                 calcite_hard_fraction = 0.02,
                 opal_hard_fraction    = 0.0,
                 dust_hard_fraction    = 0.98,

                 calcite_hard_carbon_ratio = 0.01,
                 opal_hard_carbon_ratio    = 0.01,
                 dust_hard_carbon_ratio    = 0.01,

                 particulate_organic_carbon_mass = 12.01,    # g/mol
                 calcite_mass                    = 100.09,   # g/mol
                 opal_mass                       = 60.08,    # g/mol
                 dust_mass                       = 1.0e3,    # mmol/mol — dust is measured by mass, so
                                                             # this is a unit conversion, not a molar mass

                 scale_length_depths = (100.0, 250.0, 500.0, 1000.0),   # m
                 scale_length_values = (1.0, 3.6, 4.7, 4.8),

                 oxygen_scaling_maximum = 45.0,   # mmol/m³
                 oxygen_scaling_minimum = 5.0,    # mmol/m³
                 oxygen_scaling_factor  = 2.6,

                 iron_desorption_rate     = 1.0e-4,   # per metre
                 iron_to_carbon_reference = 3.0e-6,
                 dust_to_iron             = 0.035 / 55.845 * 1.0e3,   # mmol Fe per g dust

                 surface_dust_flux     = 0.0,   # g/m²/s
                 sedimentary_iron_flux = 0.0,   # mmol/m²/s
                 open_bottom = true)

    Mc = particulate_organic_carbon_mass

    calcite_ballast_ratio = calcite_hard_carbon_ratio * calcite_mass / Mc
    opal_ballast_ratio    = opal_hard_carbon_ratio    * opal_mass    / Mc
    dust_ballast_ratio    = dust_hard_carbon_ratio    * dust_mass    / Mc

    SL = NTuple{4, FT}

    return Ballast{FT, SL, typeof(surface_dust_flux), typeof(sedimentary_iron_flux)}(
        convert(FT, particulate_organic_dissolution_length),
        convert(FT, calcite_dissolution_length),
        convert(FT, opal_dissolution_length),
        convert(FT, dust_dissolution_length),
        convert(FT, hard_dissolution_length),
        convert(FT, dust_hard_dissolution_length),
        convert(FT, calcite_hard_fraction),
        convert(FT, opal_hard_fraction),
        convert(FT, dust_hard_fraction),
        convert(FT, calcite_ballast_ratio),
        convert(FT, opal_ballast_ratio),
        convert(FT, dust_ballast_ratio),
        convert(FT, particulate_organic_carbon_mass),
        convert(FT, calcite_mass),
        convert(FT, opal_mass),
        SL(map(FT, scale_length_depths)),
        SL(map(FT, scale_length_values)),
        convert(FT, oxygen_scaling_maximum),
        convert(FT, oxygen_scaling_minimum),
        convert(FT, oxygen_scaling_factor),
        convert(FT, iron_desorption_rate),
        convert(FT, iron_to_carbon_reference),
        convert(FT, dust_to_iron),
        surface_dust_flux,
        sedimentary_iron_flux,
        open_bottom)
end

@inline lerp(x, x₀, x₁, y₀, y₁) = y₀ + (y₁ - y₀) * (x - x₀) / (x₁ - x₀)

# the depth scale-length factor: piecewise linear through the given points and flat outside. `depth` is
# positive metres below the surface, taken at the cell's LOWER FACE.
#
# Branchless: each `ifelse` overwrites the previous, so the last satisfied condition wins, and because
# `depth ≥ z[n]` implies `depth ≥ z[n-1]` that selects exactly the right segment. Every branch value is
# finite, since the depths are distinct.
@inline function scale_length(depth, b::Ballast)
    z = b.scale_length_depths
    v = b.scale_length_values

    σ = v[1]
    σ = ifelse(depth >= z[1], lerp(depth, z[1], z[2], v[1], v[2]), σ)
    σ = ifelse(depth >= z[2], lerp(depth, z[2], z[3], v[2], v[3]), σ)
    σ = ifelse(depth >= z[3], lerp(depth, z[3], z[4], v[3], v[4]), σ)
    σ = ifelse(depth >= z[4], v[4], σ)

    return σ
end

# a longer length scale, and so slower remineralisation, in low-oxygen water
@inline function oxygen_scale_factor(O₂, b::Ballast)
    hi, lo, s = b.oxygen_scaling_maximum, b.oxygen_scaling_minimum, b.oxygen_scaling_factor

    stretch = one(O₂) + (s - one(O₂)) * min(one(O₂), (hi - O₂) / (hi - lo))

    return ifelse(O₂ < hi, stretch, one(O₂))
end
