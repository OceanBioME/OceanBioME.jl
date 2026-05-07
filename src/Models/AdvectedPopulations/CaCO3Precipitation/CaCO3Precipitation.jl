"""
Abiotic calcium carbonate precipitation and dissolution model.

This module provides `SimpleCaCO3Precipitation`: a standalone biogeochemical
model that couples particulate CaCO₃ with dissolved inorganic carbon (DIC) and
total alkalinity (Alk) through kinetic precipitation under supersaturation and
dissolution under undersaturation.

Tracers
=======
* DIC   — dissolved inorganic carbon          (mmol C / m³)
* Alk   — total alkalinity                    (meq / m³)
* CaCO₃ — particulate calcium carbonate       (mmol C / m³)
* T     — temperature (°C)
* S     — salinity (PSU)

Kinetics
========
Saturation state:

    Ω = [Ca²⁺][CO₃²⁻] / Ksp

where [CO₃²⁻] is computed from DIC and Alk via carbonate equilibrium chemistry.

Precipitation rate (Ω > 1):

    R_precip = k · (Ω − 1)^n

Dissolution rate (Ω < 1):

    R_diss = k_d · (1 − Ω)^m · CaCO₃

Tracer tendencies:

    dCaCO₃/dt = +R_precip − R_diss
    dDIC/dt   = −R_precip + R_diss
    dAlk/dt   = −2 R_precip + 2 R_diss

References
==========
Morse & Arvidson (2002), Chem. Rev. — CaCO₃ dissolution kinetics
Dickson et al. (2007), PICES Special Publication 3 — carbonate chemistry
Hashim et al. (2025), Biogeosciences — CaCO₃ precipitation rate law
"""
module CaCO3PrecipitationModel

export SimpleCaCO3Precipitation, CaCO3Precipitation

using Oceananigans.Units
using Oceananigans.Fields: ZeroField, CenterField
using Oceananigans.Grids: AbstractGrid, znode, Center
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: fields
using Oceananigans.Utils: launch!

using KernelAbstractions: @kernel, @index

using OceanBioME: Biogeochemistry, ScaleNegativeTracers, setup_velocity_fields, show_sinking_velocities
using OceanBioME.Models: CarbonChemistryModel
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity,
                                     update_biogeochemical_state!

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt
import Base: show, summary

# ---------------------------------------------------------------------------
# Model struct
# ---------------------------------------------------------------------------

"""
    SimpleCaCO3Precipitation{FT, CC, CS, W}

Abiotic CaCO₃ precipitation/dissolution biogeochemical model.

Fields
======
- `precipitation_rate_constant`  — k,   mmol C m⁻³ s⁻¹
- `precipitation_order`          — n,   dimensionless
- `dissolution_rate_constant`    — k_d, s⁻¹
- `dissolution_order`            — m,   dimensionless
- `carbon_chemistry`             — `CarbonChemistry` instance for computing Ω
- `calcite_saturation`           — `CenterField` storing Ω (updated each timestep)
- `sinking_velocities`           — named tuple with sinking velocity for `:CaCO₃`
"""
struct SimpleCaCO3Precipitation{FT, CC, CS, W} <: AbstractBiogeochemistry
    precipitation_rate_constant :: FT   # k,   mmol C m⁻³ s⁻¹
    precipitation_order         :: FT   # n,   dimensionless
    dissolution_rate_constant   :: FT   # k_d, s⁻¹
    dissolution_order           :: FT   # m,   dimensionless
    carbon_chemistry            :: CC
    calcite_saturation          :: CS   # CenterField for Ω
    sinking_velocities          :: W    # m s⁻¹
end

const CaCO3Precipitation = SimpleCaCO3Precipitation

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

"""
    SimpleCaCO3Precipitation(; grid,
                               precipitation_rate_constant = 1e-5 / day,
                               precipitation_order = 2.0,
                               dissolution_rate_constant = 1e-4 / day,
                               dissolution_order = 1.0,
                               carbon_chemistry = CarbonChemistry(),
                               sinking_speed = 10 / day,
                               open_bottom = true,
                               scale_negatives = false,
                               sediment = nothing,
                               particles = nothing,
                               modifiers = nothing)

Construct a `SimpleCaCO3Precipitation` biogeochemical model with default parameters.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on
- `precipitation_rate_constant`: rate constant *k* for supersaturation-driven
  precipitation in mmol C m⁻³ s⁻¹ when (Ω − 1) = 1.  Default ≈ 10⁻⁵/day.
- `precipitation_order`: reaction order *n* for precipitation. Default 2.
- `dissolution_rate_constant`: first-order dissolution rate constant *k_d* in
  s⁻¹. Default ≈ 10⁻⁴/day.
- `dissolution_order`: reaction order *m* for dissolution. Default 1.
- `carbon_chemistry`: a `CarbonChemistry` instance used to compute Ω.
- `sinking_speed`: sinking speed of particulate CaCO₃ in m s⁻¹ (positive
  downward by convention). Default 10 m/day.
- `open_bottom`: allow particles to leave through the bottom boundary.
- `scale_negatives`: apply `ScaleNegativeTracers` modifier.
- `sediment`, `particles`, `modifiers`: optional OceanBioME slots.

Example
=======

```jldoctest
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(1, 1, 20), extent=(1, 1, 200), topology=(Periodic, Periodic, Bounded));

julia> bgc = SimpleCaCO3Precipitation(; grid)
SimpleCaCO3Precipitation{Float64} model
 Precipitation: k = 1.157e-10 mol C m⁻³ s⁻¹, n = 2.0
 Dissolution:   k_d = 1.157e-09 s⁻¹, m = 1.0
 Sinking: CaCO₃ at 1.157e-04 to 1.157e-04 m/s
```
"""
function SimpleCaCO3Precipitation(; grid::AbstractGrid{FT},
                                    precipitation_rate_constant = 1e-5 / day,
                                    precipitation_order = 2.0,
                                    dissolution_rate_constant = 1e-4 / day,
                                    dissolution_order = 1.0,
                                    carbon_chemistry = nothing,
                                    sinking_speed = 10 / day,
                                    open_bottom = true,
                                    scale_negatives = false,
                                    invalid_fill_value = NaN,
                                    sediment = nothing,
                                    particles = nothing,
                                    modifiers = nothing) where FT

    precipitation_rate_constant = convert(FT, precipitation_rate_constant)
    precipitation_order         = convert(FT, precipitation_order)
    dissolution_rate_constant   = convert(FT, dissolution_rate_constant)
    dissolution_order           = convert(FT, dissolution_order)

    carbon_chemistry = isnothing(carbon_chemistry) ? CarbonChemistry(FT) : carbon_chemistry

    sinking_speeds = (CaCO₃ = sinking_speed,)
    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    calcite_saturation_field = CenterField(grid)

    underlying_biogeochemistry =
        SimpleCaCO3Precipitation(precipitation_rate_constant,
                                 precipitation_order,
                                 dissolution_rate_constant,
                                 dissolution_order,
                                 carbon_chemistry,
                                 calcite_saturation_field,
                                 sinking_velocities)

    if scale_negatives
        scaler = ScaleNegativeTracers(underlying_biogeochemistry, grid; invalid_fill_value)
        modifiers = isnothing(modifiers) ? scaler : (modifiers..., scaler)
    end

    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation = nothing,
                           sediment,
                           particles,
                           modifiers)
end

# ---------------------------------------------------------------------------
# Biogeochemical interface
# ---------------------------------------------------------------------------

@inline required_biogeochemical_tracers(::SimpleCaCO3Precipitation) = (:DIC, :Alk, Symbol("CaCO₃"), :T, :S)

@inline required_biogeochemical_auxiliary_fields(::SimpleCaCO3Precipitation) = (:Ω,)

@inline biogeochemical_auxiliary_fields(bgc::SimpleCaCO3Precipitation) =
    (Ω = bgc.calcite_saturation,)

# ---------------------------------------------------------------------------
# Tendencies (discrete form)
# ---------------------------------------------------------------------------

@inline function precipitation_rate(bgc::SimpleCaCO3Precipitation, Ω)
    k = bgc.precipitation_rate_constant
    n = bgc.precipitation_order
    return k * max(zero(Ω), Ω - one(Ω)) ^ n
end

@inline function dissolution_rate(bgc::SimpleCaCO3Precipitation, Ω, CaCO₃)
    kd   = bgc.dissolution_rate_constant
    m    = bgc.dissolution_order
    diss = kd * max(zero(Ω), one(Ω) - Ω) ^ m * max(zero(CaCO₃), CaCO₃)
    return diss
end

# DIC: decreases by 1 mol per mol CaCO₃ precipitated, increases by 1 mol per mol dissolved
@inline function (bgc::SimpleCaCO3Precipitation)(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields)
    Ω    = @inbounds auxiliary_fields.Ω[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]

    R_precip = precipitation_rate(bgc, Ω)
    R_diss   = dissolution_rate(bgc, Ω, CaCO₃)

    return -R_precip + R_diss
end

# Alk: decreases by 2 eq per mol CaCO₃ precipitated, increases by 2 eq per mol dissolved
@inline function (bgc::SimpleCaCO3Precipitation)(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields)
    Ω    = @inbounds auxiliary_fields.Ω[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]

    R_precip = precipitation_rate(bgc, Ω)
    R_diss   = dissolution_rate(bgc, Ω, CaCO₃)

    return -2 * R_precip + 2 * R_diss
end

# CaCO₃ particulate: gains from precipitation, loses to dissolution
@inline function (bgc::SimpleCaCO3Precipitation)(i, j, k, grid, ::Val{Symbol("CaCO₃")}, clock, fields, auxiliary_fields)
    Ω     = @inbounds auxiliary_fields.Ω[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]

    R_precip = precipitation_rate(bgc, Ω)
    R_diss   = dissolution_rate(bgc, Ω, CaCO₃)

    return R_precip - R_diss
end

# Fallback: any unhandled tracer has zero BGC tendency
@inline (bgc::SimpleCaCO3Precipitation)(i, j, k, grid, ::Val, clock, fields, auxiliary_fields) =
    zero(eltype(grid))

# ---------------------------------------------------------------------------
# Saturation state update
# ---------------------------------------------------------------------------

function update_biogeochemical_state!(model, bgc::SimpleCaCO3Precipitation)
    compute_calcite_saturation!(bgc.carbon_chemistry, bgc.calcite_saturation, model)
    return nothing
end

function compute_calcite_saturation!(carbon_chemistry, calcite_saturation, model)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_calcite_saturation!,
            carbon_chemistry, calcite_saturation, grid, fields(model),
            model.clock.time)

    fill_halo_regions!(calcite_saturation, model.grid, model.clock, fields(model))

    return nothing
end

@kernel function _compute_calcite_saturation!(carbon_chemistry, calcite_saturation, grid, model_fields, t)
    i, j, k = @index(Global, NTuple)

    T   = @inbounds model_fields.T[i, j, k]
    S   = @inbounds model_fields.S[i, j, k]
    DIC = @inbounds model_fields.DIC[i, j, k]
    Alk = @inbounds model_fields.Alk[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    # Approximate hydrostatic pressure (bar): ρ₀ g |z| / 1e5
    # If we are using a box model, then set P=0, otherwise compute the pressure from the vertical coordinate.
    P = z === nothing ? zero(T) : abs(z) * 1026 * 9.80665 / 1e5

    @inbounds calcite_saturation[i, j, k] =
        CarbonChemistryModel.calcite_saturation(carbon_chemistry; DIC, T, S, Alk, P)
end

# ---------------------------------------------------------------------------
# Sinking velocity
# ---------------------------------------------------------------------------

@inline function biogeochemical_drift_velocity(bgc::SimpleCaCO3Precipitation, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end

@inline maximum_sinking_velocity(bgc::SimpleCaCO3Precipitation) =
    maximum(abs, bgc.sinking_velocities.CaCO₃)

# ---------------------------------------------------------------------------
# Adapt (GPU support)
# ---------------------------------------------------------------------------

adapt_structure(to, bgc::SimpleCaCO3Precipitation) =
    SimpleCaCO3Precipitation(adapt(to, bgc.precipitation_rate_constant),
                             adapt(to, bgc.precipitation_order),
                             adapt(to, bgc.dissolution_rate_constant),
                             adapt(to, bgc.dissolution_order),
                             adapt(to, bgc.carbon_chemistry),
                             adapt(to, bgc.calcite_saturation),
                             adapt(to, bgc.sinking_velocities))

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

summary(::SimpleCaCO3Precipitation{FT}) where FT =
    "SimpleCaCO3Precipitation{$FT} model"

function show(io::IO, bgc::SimpleCaCO3Precipitation)
    print(io, summary(bgc), "\n",
          " Precipitation: k = $(bgc.precipitation_rate_constant) mol C m⁻³ s⁻¹, ",
          "n = $(bgc.precipitation_order)\n",
          " Dissolution:   k_d = $(bgc.dissolution_rate_constant) s⁻¹, ",
          "m = $(bgc.dissolution_order)\n",
          " Sinking: CaCO₃ at ",
          show_sinking_velocities(bgc.sinking_velocities))
end

end # module
