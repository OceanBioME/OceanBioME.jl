"""
    SimpleIron(; scavenging_rate = 0.19/(360*86400),
                 ligand_stability = 1e11,
                 total_ligand = 1e-3,
                 free_iron_maximum = 3e-4)

A simplified iron cycle for the `iron` slot of [`Nutrients`](@ref). It tracks a single dissolved
iron tracer (`Fe`, in mmol Fe/m³) with two processes beyond the standard biological cycling:

1. **Ligand–iron equilibrium**: total iron is partitioned into ligand-bound and free iron via
   a stability-constant equilibrium (following Parekh et al., 2005 and MITgcm DIC `fe_chem.F`).
2. **Scavenging**: free (unbound) iron is removed at a first-order rate representing adsorption
   onto sinking particles and permanent loss to the sediment.

Free iron that exceeds `free_iron_maximum` is also removed (equivalent to MITgcm's `MINFE`
option).

External iron sources (aeolian dust deposition, sediment flux, etc.) should be applied through
Oceananigans `Forcing` on the `Fe` tracer.

Keyword Arguments
=================

- `scavenging_rate`: first-order scavenging rate of free iron (1/s)
- `ligand_stability`: ligand–free-iron stability constant (m³/mmol)
- `total_ligand`: uniform total ligand concentration (mmol/m³)
- `free_iron_maximum`: maximum concentration of free iron before precipitation (mmol/m³)
"""
@kwdef struct SimpleIron{FT}
     scavenging_rate :: FT = 0.19 / (360 * 86400)  # 1/s
    ligand_stability :: FT = 1e11                    # m³/mmol
        total_ligand :: FT = 1e-3                    # mmol/m³
   free_iron_maximum :: FT = 3e-4                    # mmol/m³
end

Adapt.adapt_structure(to, si::SimpleIron) =
    SimpleIron(adapt(to, si.scavenging_rate),
               adapt(to, si.ligand_stability),
               adapt(to, si.total_ligand),
               adapt(to, si.free_iron_maximum))

Base.summary(::SimpleIron) = string("SimpleIron (:Fe)")
function Base.show(io::IO, si::SimpleIron)
    print(io, summary(si), "\n",
          "├── Scavenging rate: ", si.scavenging_rate, "/s\n",
          "├── Ligand stability: ", si.ligand_stability, " m³/mmol\n",
          "├── Total ligand: ", si.total_ligand, " mmol/m³\n",
          "└── Free iron maximum: ", si.free_iron_maximum, " mmol/m³")
end

required_biogeochemical_tracers(::SimpleIron) = (:Fe, )

const SimpleIronNPD{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:Any, <:Any, <:SimpleIron}}

@inline function free_iron(Fe, iron::SimpleIron)
    β  = iron.ligand_stability
    Lₜ = iron.total_ligand

    discriminant = (β * Fe - β * Lₜ + 1)^2 + 4 * β * Lₜ
    L = (-β * Fe + β * Lₜ - 1 + sqrt(discriminant)) / (2 * β)

    FeL = Lₜ - L

    freefe = ifelse(Fe > 0, Fe - FeL, zero(Fe))

    return min(freefe, iron.free_iron_maximum)
end

@inline function (bgc::SimpleIronNPD)(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields)
    iron = bgc.nutrients.iron

    biological = (
        inorganic_iron_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      + inorganic_iron_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
      - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
    )

    Fe = @inbounds fields.Fe[i, j, k]
    scavenging = - iron.scavenging_rate * free_iron(Fe, iron)

    return biological + scavenging
end
