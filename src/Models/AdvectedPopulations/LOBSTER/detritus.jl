using Oceananigans.Fields: ZeroField, ConstantField

"""
    TwoParticleAndDissolved

`TwoParticleAndDissolved` defines the default living component for the `LOBSTER` 
biogeochemical model. It includes small and large particulate organic matter 
(`sPOM` and `bPOM`), and dissolved organic matter (`DOM`).

`POM` is produced by "waste" from the biological system (e.g. excretion from
zooplankton or phytoplankton death) and is portioned into `sPOM` and `bPOM` 
as by the `small_solid_waste_fraction` parameter. `sPOM` is depleated by grazing
by the `biology` (typically the zooplankton), and is remineralised. Calcite 
produced by the `biology` (from phytoplankton mortality) is assumed to all end up
in `bPOM`. 

`DOM` is produced by the biological system (e.g. excretion) and the breakdown of 
`POM`, and is broken down at a constant remineralisation rate. 

Note: we refer to the particle classes as "small" and "large", but label them 
`sPOM` and `bPOM` for "small" and "big" as using `lPOM` may be confused for
labile `POM`. We agree that this is a clumsy solution and would be happy to find
an alternative.
""" 
@kwdef struct TwoParticleAndDissolved{FT, SS, LS}
    remineralisation_inorganic_fraction :: FT = 0.0     # 
            small_reminerslisation_rate :: FT = 5.88e-7 # 1/s
             large_reminerlisation_rate :: FT = 5.88e-7 # 1/s
         dissolved_reminerlisation_rate :: FT = 3.86e-7 # 1/s

             small_solid_waste_fraction :: FT = 0.5     # 

                         redfield_ratio :: FT = 6.56    # mol C/mol N

        small_particle_sinking_velocity :: SS = (u = ZeroField(), v = ZeroField(), w = ConstantField(-3.47e-5))
        large_particle_sinking_velocity :: LS = (u = ZeroField(), v = ZeroField(), w = ConstantField(-200/day))
end

required_biogeochemical_tracers(::TwoParticleAndDissolved) = (:sPOM, :bPOM, :DOM)

@inline small_particulate_concentration(::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.sPOM[i, j, k]

@inline large_particulate_concentration(::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.bPOM[i, j, k]

@inline dissolved_organic_nitrogen(::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) = 
    @inbounds fields.DOM[i, j, k]

@inline small_particulate_carbon_concentration(detritus::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.sPOM[i, j, k] * detritus.redfield_ratio

@inline large_particulate_carbon_concentration(detritus::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.bPOM[i, j, k] * detritus.redfield_ratio

@inline dissolved_organic_carbon(detritus::TwoParticleAndDissolved, i, j, k, fields, auxiliary_fields) = 
    @inbounds fields.DOM[i, k, k] * detritus.redfield_ratio

"""
    VariableRedfieldDetritus

`VariableRedfieldDetritus` defines builds on the default `TwoParticleAndDissolved`
component for the `LOBSTER` biogeochemical model, splitting each component into 
separate carbon and nitrogen compartements: `sPOC`, `bPOC`, `DOC`, `sPON`, `bPON`,
and `DON`. 

Standing alone this will produce identical behavior to the `TwoParticleAndDissolved`
formulation, but if particles are injected with a different C:N ratio (for example)
from the `SugarKelp` model, then the ratio will evolve differently.
""" 
@kwdef struct VariableRedfieldDetritus{FT, SS, LS}
    remineralisation_inorganic_fraction :: FT = 0.0     # 
            small_reminerslisation_rate :: FT = 5.88e-7 # 1/s
             large_reminerlisation_rate :: FT = 5.88e-7 # 1/s
         dissolved_reminerlisation_rate :: FT = 3.86e-7 # 1/s

             small_solid_waste_fraction :: FT = 0.5     # 

        small_particle_sinking_velocity :: SS = (u = ZeroField(), v = ZeroField(), w = ConstantField(-3.47e-5))
        large_particle_sinking_velocity :: LS = (u = ZeroField(), v = ZeroField(), w = ConstantField(-200/day))
end

required_biogeochemical_tracers(::VariableRedfieldDetritus) = (:sPOC, :bPOC, :DOC, :sPON, :bPON, :DON)

@inline (lobster::LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus})(i, j, k, grid, val_name::Val{:sPOC}, clock, fields, auxiliary_fields) = (
    lobster.detritus.small_solid_waste_fraction * solid_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  - grazing(lobster, i, j, k, val_name, fields, auxiliary_fields)
  - lobster.detritus.small_reminerslisation_rate * small_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

@inline (lobster::LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus})(i, j, k, grid, ::Val{:bPOC}, clock, fields, auxiliary_fields) = (
    (1 - lobster.detritus.small_solid_waste_fraction) * solid_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  + calcite_production(lobster, i, j, k, fields, auxiliary_fields)
  - lobster.detritus.large_reminerlisation_rate * large_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

@inline (lobster::LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus})(i, j, k, grid, ::Val{:DOC}, clock, fields, auxiliary_fields) = (
    biology_organic_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  + detritus_organic_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  - lobster.detritus.dissolved_reminerlisation_rate * dissolved_organic_carbon(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

@inline small_particulate_concentration(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.sPON[i, j, k]

@inline large_particulate_concentration(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.bPON[i, j, k]

@inline dissolved_organic_nitrogen(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) = 
    @inbounds fields.DON[i, k, k]

@inline small_particulate_carbon_concentration(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.sPOC[i, j, k]

@inline large_particulate_carbon_concentration(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.bPOC[i, j, k]

@inline dissolved_organic_carbon(::VariableRedfieldDetritus, i, j, k, fields, auxiliary_fields) = 
    @inbounds fields.DOC[i, k, k]

##### common

const TWO_PARTICLE_SIZE_LOBSTER = Union{LOBSTER{<:Any, <:Any, <:TwoParticleAndDissolved}, LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus}}
const TWO_PARTICLE_SIZES = Union{TwoParticleAndDissolved, VariableRedfieldDetritus}

const SMALL_SINKING_OFF =  Union{LOBSTER{<:Any, <:Any, <:TwoParticleAndDissolved{<:Any, Nothing}}, 
                                 LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus{<:Any, Nothing}}}

const LARGE_SINKING_OFF =  Union{LOBSTER{<:Any, <:Any, <:TwoParticleAndDissolved{<:Any, <:Any, Nothing}}, 
                                 LOBSTER{<:Any, <:Any, <:VariableRedfieldDetritus{<:Any, <:Any, Nothing}}}

const SMALL_PARTICLES = Union{Val{:sPOM}, Val{:sPON}, Val{:sPOC}}
const LARGE_PARTICLES = Union{Val{:bPOM}, Val{:bPON}, Val{:bPOC}}

biogeochemical_drift_velocity(lobster::TWO_PARTICLE_SIZE_LOBSTER, ::SMALL_PARTICLES) = 
    lobster.detritus.small_particle_sinking_velocity
biogeochemical_drift_velocity(lobster::TWO_PARTICLE_SIZE_LOBSTER, ::LARGE_PARTICLES) = 
    lobster.detritus.large_particle_sinking_velocity

biogeochemical_drift_velocity(::SMALL_SINKING_OFF, ::SMALL_PARTICLES) = nothing
biogeochemical_drift_velocity(::LARGE_SINKING_OFF, ::LARGE_PARTICLES) = nothing

@inline (lobster::TWO_PARTICLE_SIZE_LOBSTER)(i, j, k, grid, val_name::Union{Val{:DOM}, Val{:DON}}, clock, fields, auxiliary_fields) = (
    biology_organic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  + detritus_organic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  - lobster.detritus.dissolved_reminerlisation_rate * dissolved_organic_nitrogen(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

@inline (lobster::TWO_PARTICLE_SIZE_LOBSTER)(i, j, k, grid, val_name::Union{Val{:sPOM}, Val{:sPON}}, clock, fields, auxiliary_fields) = (
    lobster.detritus.small_solid_waste_fraction * solid_waste(lobster, i, j, k, fields, auxiliary_fields)
  - grazing(lobster, i, j, k, val_name, fields, auxiliary_fields)
  - lobster.detritus.small_reminerslisation_rate * small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

@inline (lobster::TWO_PARTICLE_SIZE_LOBSTER)(i, j, k, grid, val_name::Union{Val{:bPOM}, Val{:bPON}}, clock, fields, auxiliary_fields) = (
    (1 - lobster.detritus.small_solid_waste_fraction) * solid_waste(lobster, i, j, k, fields, auxiliary_fields)
  - lobster.detritus.large_reminerlisation_rate * large_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
)

##### waste
@inline function detritus_inorganic_nitrogen_waste(lobster::TWO_PARTICLE_SIZE_LOBSTER, i, j, k, fields, auxiliary_fields)
    α  = lobster.detritus.remineralisation_inorganic_fraction
    sμ = lobster.detritus.small_reminerslisation_rate
    bμ = lobster.detritus.large_reminerlisation_rate
    dμ = lobster.detritus.dissolved_reminerlisation_rate

    sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    bPOM = large_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    DOM  = dissolved_organic_nitrogen(lobster.detritus, i, j, k, fields, auxiliary_fields)

    return (α * (sμ * sPOM + bμ * bPOM) + dμ * DOM)
end

@inline function detritus_organic_nitrogen_waste(lobster::TWO_PARTICLE_SIZE_LOBSTER, i, j, k, fields, auxiliary_fields)
    α  = lobster.detritus.remineralisation_inorganic_fraction
    sμ = lobster.detritus.small_reminerslisation_rate
    bμ = lobster.detritus.large_reminerlisation_rate

    sPOM = small_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    bPOM = large_particulate_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)

    return (1 - α) * (sμ * sPOM + bμ * bPOM)
end

@inline function detritus_inorganic_carbon_waste(lobster::TWO_PARTICLE_SIZE_LOBSTER, i, j, k, fields, auxiliary_fields)
    α  = lobster.detritus.remineralisation_inorganic_fraction
    sμ = lobster.detritus.small_reminerslisation_rate
    bμ = lobster.detritus.large_reminerlisation_rate
    dμ = lobster.detritus.dissolved_reminerlisation_rate

    sPOC = small_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    bPOC = large_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    DOC  = dissolved_organic_carbon(lobster.detritus, i, j, k, fields, auxiliary_fields)

    return α * (sμ * sPOC + bμ * bPOC) + dμ * DOC
end

@inline function detritus_organic_carbon_waste(lobster::TWO_PARTICLE_SIZE_LOBSTER, i, j, k, fields, auxiliary_fields)
    α  = lobster.detritus.remineralisation_inorganic_fraction
    sμ = lobster.detritus.small_reminerslisation_rate
    bμ = lobster.detritus.large_reminerlisation_rate

    sPOC = small_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)
    bPOC = large_particulate_carbon_concentration(lobster.detritus, i, j, k, fields, auxiliary_fields)

    return (1 - α) * (sμ * sPOC + bμ * bPOC)
end

"""
    TwoParticleAndDissolved(grid; 
                            small_particle_sinking_speed = 3.47e-5, # m/s
                            large_particle_sinking_speed = 200/day, # m/s
                            open_bottom = true,
                            kwargs...)

Construct an instance of `TwoParticleAndDissolved` specifying the 
`small_particle_sinking_speed` and `large_particle_sinking_speed`. If `open_bottom`
is true then particles sink out of the bottom, if false then the sinking velocity 
smoothly goes to zero at the bottom to prevent the tracers leaving the domain.
"""
function TwoParticleAndDissolved(grid; 
                                 small_particle_sinking_speed = 3.47e-5, # m/s
                                 large_particle_sinking_speed = 200/day, # m/s
                                 open_bottom = true,
                                 kwargs...)


    small_particle_sinking_velocity = setup_velocity_fields((; sPOM = small_particle_sinking_speed), grid, open_bottom; three_D = true).sPOM
    large_particle_sinking_velocity = setup_velocity_fields((; bPOM = large_particle_sinking_speed), grid, open_bottom; three_D = true).bPOM

    return TwoParticleAndDissolved(;small_particle_sinking_velocity, large_particle_sinking_velocity, kwargs...)
end

"""
    VariableRedfieldDetritus(grid; 
                             small_particle_sinking_speed = 3.47e-5, # m/s
                             large_particle_sinking_speed = 200/day, # m/s
                             open_bottom = true,
                             kwargs...)

Construct an instance of `TwoParticleAndDissolved` specifying the 
`small_particle_sinking_speed` and `large_particle_sinking_speed`. If `open_bottom`
is true then particles sink out of the bottom, if false then the sinking velocity 
smoothly goes to zero at the bottom to prevent the tracers leaving the domain.
"""
function VariableRedfieldDetritus(grid; 
                                  small_particle_sinking_speed = 3.47e-5, # m/s
                                  large_particle_sinking_speed = 200/day, # m/s
                                  open_bottom = true,
                                  kwargs...)


    small_particle_sinking_velocity = setup_velocity_fields((; sPOM = small_particle_sinking_speed), grid, open_bottom; three_D = true).sPOM
    large_particle_sinking_velocity = setup_velocity_fields((; bPOM = large_particle_sinking_speed), grid, open_bottom; three_D = true).bPOM

    return VariableRedfieldDetritus(;small_particle_sinking_velocity, large_particle_sinking_velocity, kwargs...)
end