# the default assumptions

function dissolved_waste end
function solid_waste end
function calcium_carbonate_dissolution end
function inorganic_waste end
function nutrient_uptake end

@inline nitrogen_ratio(i, j, k, grid, plankton, bgc, fields) = nitrogen_ratio(plankton, bgc)
@inline nitrogen_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = one(FT)

@inline carbon_ratio(i, j, k, grid, plankton, bgc, fields) = carbon_ratio(plankton, bgc)
@inline carbon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 106/16)

@inline phosphate_ratio(i, j, k, grid, plankton, bgc, fields) = phosphate_ratio(plankton, bgc)
@inline phosphate_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 1/16)

@inline iron_ratio(i, j, k, grid, plankton, bgc, fields) = iron_ratio(plankton, bgc)
@inline iron_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 0.0032/16) 

@inline silicon_ratio(i, j, k, grid, plankton, bgc, fields) = silicon_ratio(plankton, bgc)
@inline silicon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = zero(FT)

@inline calcium_carbonate_rain_ratio(i, j, k, grid, plankton, bgc, fields) = calcium_carbonate_rain_ratio(plankton, bgc)
@inline calcium_carbonate_rain_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = zero(FT)

@inline chlorophyll_ratio(plankton, ::Val) = chlorophyll_ratio(plankton)

# Organic matter carried down by the inorganic carbon's particles (MARBL's "hard" POC, ballasted by the CaCO₃), in
# the plankton's units: its production is taken out of the particulate detritus production and its
# remineralisation is returned with the detritus' (see `DissolvedParticulate`). Only `ImplicitCalciumCarbonate`
# carries any.
@inline hard_organic_production(i, j, k, grid, inorganic_carbon, ::NutrientsPlanktonDetritus{FT}, fields, auxiliary_fields) where FT = zero(FT)
@inline hard_organic_remineralisation(i, j, k, grid, inorganic_carbon, ::NutrientsPlanktonDetritus{FT}, fields, auxiliary_fields) where FT = zero(FT)

# The sinking particle mass flux into the top of cell `k` (ng cm⁻² s⁻¹, MARBL's units, positive downward), which
# `SimpleIron`'s `ParticleScavenging` scales the scavenging with. Only components that sink implicitly and store their
# flux (`store_flux = true`) contribute; `stores_sinking_flux` says whether one does.
@inline sinking_particle_mass_flux(i, j, k, grid, component, ::NutrientsPlanktonDetritus{FT}) where FT = zero(FT)
stores_sinking_flux(component) = false

# MARBL's mass of sinking organic matter, 3 × its carbon (g / mol C), and of calcium carbonate (g / mol), and the
# conversion of a flux from mmol m⁻² s⁻¹ to nmol cm⁻² s⁻¹ (so that mol × g mol⁻¹ gives ng cm⁻² s⁻¹)
const organic_matter_mass_per_carbon = 3 * 12.01
const calcium_carbonate_molar_mass   = 100.09
const mmol_m⁻²_to_nmol_cm⁻²          = 100

# the total of the stored implicit sinking fluxes (`F[i, j, k]` leaves the bottom of cell `k`) into the top of cell `k`
@inline sinking_flux_into(i, j, k, grid, flux::NamedTuple) =
    sum(map(F -> ifelse(k == size(grid, 3), zero(eltype(grid)), @inbounds(F[i, j, k + 1])), values(flux)))

for (element, symbol) in pairs((nitrogen = :N, phosphate = :PO₄, iron = :Fe, silicon = :Si, carbon = :C))
    inorganic_waste_name  = Symbol(:inorganic_, element, :_waste)
    solid_waste_name      = Symbol(:solid_,     element, :_waste)
    dissolved_waste_name  = Symbol(:dissolved_, element, :_waste)
    ratio_name            = Symbol(element, :_ratio)
    @eval begin
        @inline $inorganic_waste_name(i, j, k, grid, plankton_or_detritus, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            inorganic_waste(i, j, k, grid, plankton_or_detritus, bgc, fields, auxiliary_fields)

        @inline $solid_waste_name(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            solid_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

        @inline $dissolved_waste_name(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            dissolved_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

        @inline nutrient_uptake(i, j, k, grid, ::Val{$(QuoteNode(symbol))}, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
    end
end
