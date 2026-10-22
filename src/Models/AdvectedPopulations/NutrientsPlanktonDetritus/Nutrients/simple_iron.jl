"""
    SimpleIron(; scavenging_rate = 0.19/(360*86400),
                 scavenging = scavenging_rate,
                 ligand_stability = 1e5,
                 total_ligand = 1e-3)

A simplified iron cycle for the `iron` slot of [`Nutrients`](@ref). It tracks a single dissolved
iron tracer (`Fe`, in mmol Fe/m³) with two processes beyond the standard biological cycling:

1. **Ligand–iron equilibrium**: total iron is partitioned into ligand-bound and free iron via
   a stability-constant equilibrium (following Parekh et al., 2005 and MITgcm DIC `fe_chem.F`).
2. **Scavenging**: iron adsorbed onto sinking particles and permanently lost. By default (a number,
   `scavenging_rate`) free iron is removed at that first-order rate, as MITgcm DIC does; with
   [`ParticleScavenging`](@ref) free and ligand-bound iron are removed at rates proportional to the sinking
   particle mass, as MARBL does.

Free iron is not capped (as MITgcm DIC without its `MINFE` option), so iron above the total ligand
is scavenged in proportion to the excess.

External iron sources (aeolian dust deposition, sediment flux, etc.) should be applied through
Oceananigans `Forcing` on the `Fe` tracer.

Keyword Arguments
=================

- `scavenging_rate`: first-order scavenging rate of free iron (1/s)
- `scavenging`: the scavenging, either that rate or a [`ParticleScavenging`](@ref)
- `ligand_stability`: ligand–free-iron stability constant (m³/mmol; MITgcm's `ligand_stab = 1e8` m³/mol)
- `total_ligand`: uniform total ligand concentration (mmol/m³)
"""
struct SimpleIron{FT, S}
          scavenging :: S
    ligand_stability :: FT
        total_ligand :: FT
end

"""
    ParticleScavenging(FT = Float64; free_rate = 22, bound_rate = 1.2)

MARBL's iron scavenging for [`SimpleIron`](@ref) (`compute_scavenging` in `marbl_interior_tendency_mod.F90`): free
and ligand-bound iron are removed at

    (free_rate × Fe′ + bound_rate × FeL) × M  per year,

where `M` is the mass flux of sinking particles into the top of the cell in ng cm⁻² s⁻¹, so that the rates are per
year per ng cm⁻² s⁻¹ (MARBL's `parm_Fe_scavenge_rate0` and `parm_FeLig_scavenge_rate0`). As in MARBL, `M` counts the
organic matter as 3 × its carbon mass and calcium carbonate at 100.09 g mol⁻¹ with the organic matter it carries,
taken from the stored implicit sinking fluxes of the detritus and inorganic carbon (`store_flux = true`): the
detritus must store its flux, and inorganic carbon without sinking calcium carbonate adds nothing. Unlike MARBL,
biogenic silica and dust are not counted, and the scavenged iron is lost rather than carried down by the particles.
"""
struct ParticleScavenging{FT}
     free_rate :: FT # yr⁻¹ (ng cm⁻² s⁻¹)⁻¹
    bound_rate :: FT # yr⁻¹ (ng cm⁻² s⁻¹)⁻¹
end

ParticleScavenging(FT = Float64; free_rate = 22, bound_rate = 1.2) =
    ParticleScavenging(convert(FT, free_rate), convert(FT, bound_rate))

convert_scavenging(FT, rate::Number) = convert(FT, rate)
convert_scavenging(FT, s::ParticleScavenging) = ParticleScavenging(FT; free_rate = s.free_rate, bound_rate = s.bound_rate)

SimpleIron{FT}(; scavenging_rate = 0.19 / (360 * 86400), # 1/s
                 scavenging = scavenging_rate,
                 ligand_stability = 1e5,                 # m³/mmol
                 total_ligand = 1e-3) where FT =         # mmol/m³
    SimpleIron(convert_scavenging(FT, scavenging), convert(FT, ligand_stability), convert(FT, total_ligand))

SimpleIron(; kwargs...) = SimpleIron{Float64}(; kwargs...)

Adapt.adapt_structure(to, si::SimpleIron) =
    SimpleIron(adapt(to, si.scavenging),
               adapt(to, si.ligand_stability),
               adapt(to, si.total_ligand))

Base.summary(::SimpleIron) = string("SimpleIron (:Fe)")

scavenging_summary(rate::Number) = string("free iron at ", rate, "/s")
scavenging_summary(s::ParticleScavenging) =
    string("MARBL-like, free iron ", s.free_rate, " and bound iron ", s.bound_rate, " yr⁻¹ (ng cm⁻² s⁻¹)⁻¹ of sinking particles")

function Base.show(io::IO, si::SimpleIron)
    print(io, summary(si), "\n",
          "├── Scavenging: ", scavenging_summary(si.scavenging), "\n",
          "├── Ligand stability: ", si.ligand_stability, " m³/mmol\n",
          "└── Total ligand: ", si.total_ligand, " mmol/m³")
end

required_biogeochemical_tracers(::SimpleIron) = (:Fe, )

const SimpleIronNPD{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:Any, <:Any, <:SimpleIron}}

# the free (Fe′) and ligand-bound (FeL) iron, both zero where there is no iron
@inline function free_and_bound_iron(Fe, iron::SimpleIron)
    β  = iron.ligand_stability
    Lₜ = iron.total_ligand

    discriminant = (β * Fe - β * Lₜ + 1)^2 + 4 * β * Lₜ
    L = (-β * Fe + β * Lₜ - 1 + sqrt(discriminant)) / (2 * β)

    FeL = Lₜ - L

    return ifelse(Fe > 0, Fe - FeL, zero(Fe)), ifelse(Fe > 0, FeL, zero(Fe))
end

@inline free_iron(Fe, iron::SimpleIron) = first(free_and_bound_iron(Fe, iron))

@inline iron_scavenging(i, j, k, grid, rate::Number, Fe, bgc) = rate * free_iron(Fe, bgc.nutrients.iron)

@inline function iron_scavenging(i, j, k, grid, s::ParticleScavenging, Fe, bgc::NutrientsPlanktonDetritus{FT}) where FT
    M = sinking_particle_mass_flux(i, j, k, grid, bgc.detritus, bgc) +
        sinking_particle_mass_flux(i, j, k, grid, bgc.inorganic_carbon, bgc) # ng cm⁻² s⁻¹

    Fe′, FeL = free_and_bound_iron(Fe, bgc.nutrients.iron)

    return M * (s.free_rate * Fe′ + s.bound_rate * FeL) / convert(FT, 365 * 86400) # MARBL's year
end

# `ParticleScavenging` needs the detritus' sinking flux
validate_iron_scavenging(iron, detritus) = nothing

validate_iron_scavenging(::SimpleIron{<:Any, <:ParticleScavenging}, detritus) =
    stores_sinking_flux(detritus) ||
        throw(ArgumentError("`ParticleScavenging` needs detritus that sinks implicitly and stores its flux " *
                            "(e.g. `DissolvedParticulate(grid; implicit_sinking = true, store_flux = true)`), got $(summary(detritus))"))

@inline function (bgc::SimpleIronNPD)(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields)
    iron = bgc.nutrients.iron

    biological = (
        inorganic_iron_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      + inorganic_iron_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
      - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
    )

    Fe = @inbounds fields.Fe[i, j, k]

    return biological - iron_scavenging(i, j, k, grid, iron.scavenging, Fe, bgc)
end
