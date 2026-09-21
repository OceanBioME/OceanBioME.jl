# MARBL's dissolved/refractory detritus model with implicit particle sinking tracking carbon/nitrogen/phosphate

# defaults for the two generics declared in the module files: without an oxygen component there is no
# low-oxygen stretch, and without an iron cycle nothing scavenges so the sinking mass is unused
@inline ballast_oxygen_scale(i, j, k, grid, oxygen, b, fields) = one(eltype(grid))

#####
##### Implicit ballast sinking — the column sweep
#####
##### The particulates are NOT tracers here. Their flux profile is rebuilt from scratch every step and
##### thrown away, and only the resulting remineralisation reaches the tendency equations: no storage,
##### no time lag, no sinking CFL condition.
#####
##### The governing statement is the steady-state budget
#####
#####     0 = -∂F/∂z + Π - R      ⟺      R = Π - ∂F/∂z
#####
##### with F ≤ 0 downward. Everything produced in a cell is either remineralised there or carried away
##### within the same step.
#####
##### Two discretisation choices are reproduced deliberately, because they change the numbers and the
##### parameter values are tuned around them: the depth factor is evaluated at the cell's LOWER FACE
##### rather than its centre, and the hard sub-class injects its production undecayed.
#####
struct MultiElementRefractoryDissolved{FT, R, F, B, RM, SM, FL, FP, FI} <: AbstractMultiElementRefractoryDissolvedDetritus
    semilabile_remineralisation_rate_light :: R
     semilabile_remineralisation_rate_dark :: R
          refractory_remineralisation_rate :: R

                     photodegradation_rate :: FT
                        uv_reference_depth :: FT

            refractory_production_fraction :: F
           particulate_refractory_fraction :: F

                                   ballast :: B

                          remineralisation :: RM   # NamedTuple of centre fields
                              sinking_mass :: SM   # kept because ligand scavenging does not cancel out of the ligand tendency
                                floor_flux :: FL   # positive-down magnitudes, for the sediment
                                    fluxes :: FP   # optional per-face profiles, for diagnostics only

                             floor_indices :: FI
end

"""
    MultiElementRefractoryDissolved(grid; ballast = Ballast(FT), flux_diagnostics = false)

A multi-element detritus which sinks its particulates through the implicit mineral-protection (ballast)
flux of Armstrong et al. (2000) rather than carrying them as tracers — the sibling of
[`MultiElementRefractoryDissolvedParticulate`](@ref), swappable in the `detritus` slot.

The dissolved organic matter (`DOC`/`DON`/`DOP` and refractory `DOCr`/`DONr`/`DOPr`) is identical to its
sibling. The difference is the particulates: organic carbon and phosphorus, biogenic silica, particulate
iron, dust and calcite are **not tracers**. Each step the column is swept from the surface down solving
the sinking-flux equation, and the resulting remineralisation is what the tendencies read.

Requires a [`ImplicitExplicitCalcite`](@ref) in the `inorganic_carbon` slot, since calcite is one of the three
ballast minerals and the sweep cannot run without it.

Keyword Arguments
=================

- `ballast`: the [`Ballast`](@ref) sinking parameters
- `uv_reference_depth`: the depth the surface photodegradation dose is spread over (m)
- `flux_diagnostics`: allocate per-face flux profiles for diagnostics (nothing in the model reads them)
"""
function MultiElementRefractoryDissolved(grid::AbstractGrid{FT};
                                         photodegradation_rate = 1/(18yr),
                                         uv_reference_depth = 10.0,
                                         ballast = Ballast(FT),
                                         light_rates = (carbon = 1/(15yr),    nitrogen = 1/(15yr),   phosphorus = 1/(60yr)),
                                         dark_rates = (carbon = 1/(6yr),     nitrogen = 1/(5.5yr),  phosphorus = 1/(4.5yr)),
                                         reference_rates = (carbon = 1/(16000yr), nitrogen = 1/(9500yr), phosphorus = 1/(5500yr)),
                                         production_refractory_fraction = (carbon = 0.01,      nitrogen = 0.0115,      phosphorus = 0.003),
                                         particulate_refractory_fraction = (carbon = 0.01*0.06, nitrogen = 0.0115*0.03, phosphorus = 0.003*0.06),
                                         flux_diagnostics = false) where FT

    particulates = (:POC, :POP, :bSi, :PFe) # not generalising this yet

    remineralisation = NamedTuple{particulates}(ntuple(_ -> CenterField(grid), length(particulates)))
    floor_flux       = NamedTuple{particulates}(ntuple(_ -> Field{Center, Center, Nothing}(grid),
                                                       length(particulates)))
    sinking_mass = CenterField(grid)

    flux_names = (:POC_soft, :POC_hard, :CaCO₃, :bSi, :dust, :PFe, :POP)
    fluxes = flux_diagnostics ?
        NamedTuple{flux_names}(ntuple(_ -> ZFaceField(grid), length(flux_names))) : nothing

    RT = NamedTuple{(:carbon, :nitrogen, :phosphorus), NTuple{3, FT}}

    floor_indices = floor_index_field(grid)

    return MultiElementRefractoryDissolved{FT, RT, RT, typeof(ballast),
                                          typeof(remineralisation), typeof(sinking_mass),
                                          typeof(floor_flux), typeof(fluxes), typeof(floor_indices)}(
        RT(light_rates),
        RT(dark_rates),
        RT(reference_rates),
        convert(FT, photodegradation_rate),
        convert(FT, uv_reference_depth),
        RT(production_refractory_fraction),
        RT(particulate_refractory_fraction),
        ballast,
        remineralisation,
        sinking_mass,
        floor_flux,
        fluxes,
        floor_indices
    )
end

const MultiElementMultiElementRefractoryDissolvedDetritus_NPD =
    NutrientsPlanktonDetritus{<:Any, <:Nutrients, <:Any, <:MultiElementRefractoryDissolved, <:Any}

required_biogeochemical_tracers(::MultiElementRefractoryDissolved) = 
    (:DOC, :DON, :DOP, :DOCr, :DONr, :DOPr)

required_biogeochemical_auxiliary_fields(::MultiElementRefractoryDissolved) =
    (:POC_remin, :POP_remin, :bSi_remin, :PFe_remin, :sinking_mass,
     :POC_floor, :POP_floor, :bSi_floor, :PFe_floor)

biogeochemical_auxiliary_fields(d::MultiElementRefractoryDissolved) =
    (POC_remin = d.remineralisation.POC, 
     POP_remin = d.remineralisation.POP,
     bSi_remin = d.remineralisation.bSi, 
     PFe_remin = d.remineralisation.PFe,
     sinking_mass = d.sinking_mass,
     POC_floor = d.floor_flux.POC, POP_floor = d.floor_flux.POP,
     bSi_floor = d.floor_flux.bSi, PFe_floor = d.floor_flux.PFe)

@inline poc_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolved, fields) =
    @inbounds d.remineralisation.POC[i, j, k]

@inline pop_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolved, fields) =
    @inbounds d.remineralisation.POP[i, j, k]

@inline bsi_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolved, fields) =
    @inbounds d.remineralisation.bSi[i, j, k]

@inline pfe_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolved, fields) =
    @inbounds d.remineralisation.PFe[i, j, k]

@inline scavenging_sinking_mass(i, j, k, grid, bgc::MultiElementMultiElementRefractoryDissolvedDetritus_NPD, fields, ci::ComplexedIron) =
    @inbounds bgc.detritus.sinking_mass[i, j, k]

@inline get_surface_dust_flux(i, j, grid, Φ::Number) = Φ
@inline get_surface_dust_flux(i, j, grid, Φ) = @inbounds Φ[i, j, 1]

@inline get_sedimentary_iron_flux(i, j, k, grid, Φ::Number) = Φ
@inline get_sedimentary_iron_flux(i, j, k, grid, Φ) = @inbounds Φ[i, j, k]

@inline store_fluxes!(::Nothing, i, j, k, FsPOC, FhPOC, FCaCO₃, FbSi, Fdust, FPFe, FPOP) = nothing

@inline function store_fluxes!(F, i, j, k, FsPOC, FhPOC, FCaCO₃, FbSi, Fdust, FPFe, FPOP)
    @inbounds begin
        F.POC_soft[i, j, k] = FsPOC
        F.POC_hard[i, j, k] = FhPOC
           F.CaCO₃[i, j, k] = FCaCO₃
             F.bSi[i, j, k] = FbSi
            F.dust[i, j, k] = Fdust
             F.PFe[i, j, k] = FPFe
             F.POP[i, j, k] = FPOP
    end

    return nothing
end

# the index the column sweep starts from. It has to be an integer field: the sweep uses it as a loop bound
# and compares it to `k`, so a float would make the range non-integer.
floor_index_field(grid) = OneField(Int)
function floor_index_field(grid::ImmersedBoundaryGrid)
    floor_indices = Field{Center, Center, Nothing}(grid, Int)

    launch!(architecture(grid), grid, :xy, compute_floor_indices!, grid, floor_indices, size(grid, 3))

    return floor_indices
end

@kernel function compute_floor_indices!(grid, floor_indices, Nz)
    i, j = @index(Global, NTuple)

    k = 1

    while immersed_cell(i, j, k, grid.underlying_grid, grid.immersed_boundary) && k < Nz
        k += 1
    end

    @inbounds floor_indices[i, j, 1] = k
end

@kernel function compute_ballast_fluxes!(grid, bgc, model_fields, aux)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    
    detritus  = bgc.detritus
    ballast  = detritus.ballast
    inorganic_carbon = bgc.inorganic_carbon

    k_bottom = detritus.floor_indices[i, j, 1]

    ρ_calcite = ballast.calcite_ballast_ratio
    ρ_opal = ballast.opal_ballast_ratio
    ρ_dust = ballast.dust_ballast_ratio

    γ_calcite = ballast.calcite_hard_fraction
    γ_opal = ballast.opal_hard_fraction

    ∅ = zero(grid)

    Φdust = get_surface_dust_flux(i, j, grid, ballast.surface_dust_flux)

    FsCaCO₃ = ∅; FhCaCO₃ = ∅
    FsbSi   = ∅; FhbSi   = ∅
    FsPOC   = ∅; FhPOC   = ∅
    FPOP    = ∅
    FPFe    = ∅

    Fsdust = -(one(Φdust) - ballast.dust_hard_fraction) * Φdust
    Fhdust = -ballast.dust_hard_fraction * Φdust

    𝒬 = ρ_dust * Φdust

    store_fluxes!(detritus.fluxes, i, j, Nz + 1, FsPOC, FhPOC, FsCaCO₃ + FhCaCO₃, FsbSi + FhbSi,
                  Fsdust + Fhdust, FPFe, FPOP)

    for k in Nz:-1:k_bottom
        Δz  = Δzᶜᶜᶜ(i, j, k, grid)
        Δz⁻ = one(Δz) / Δz

        depth = znode(i, j, k, grid, Center(), Center(), Face())
        σ     = scale_length(depth, ballast)
        sO₂   = ballast_oxygen_scale(i, j, k, grid, bgc.oxygen, ballast, model_fields)
        scale = σ * sO₂

        ℓPOC   = scale * ballast.particulate_organic_dissolution_length
        ℓCaCO₃ = scale * ballast.calcite_dissolution_length
        ℓbSi   = scale * ballast.opal_dissolution_length
        ℓdust  = scale * ballast.dust_dissolution_length

        dPOC   = exp(-Δz / ℓPOC)
        dCaCO₃ = exp(-Δz / ℓCaCO₃)
        dbSi   = exp(-Δz / ℓbSi)
        ddust  = exp(-Δz / ℓdust)
        dhard  = exp(-Δz / ballast.hard_dissolution_length)
        dhdust = exp(-Δz / ballast.dust_hard_dissolution_length)

        ΠPOC   = poc_production(i, j, k, grid, bgc, model_fields, aux)
        ΠPOP   = particulate_phosphorus_production(i, j, k, grid, bgc.plankton, bgc, model_fields, aux)
        ΠCaCO₃ = particulate_calcium_carbonate_production(i, j, k, grid, bgc.plankton, bgc, model_fields, aux)
        ΠbSi   = biogenic_silica_production(i, j, k, grid, bgc.plankton, bgc, model_fields, aux)

        FsCaCO₃⁺, FhCaCO₃⁺ = FsCaCO₃, FhCaCO₃
        FsbSi⁺,   FhbSi⁺   = FsbSi,   FhbSi
        Fsdust⁺,  Fhdust⁺  = Fsdust,  Fhdust
        FsPOC⁺,   FhPOC⁺   = FsPOC,   FhPOC
        FPOP⁺,    FPFe⁺    = FPOP,    FPFe

        𝓜 = ballast_sinking_mass(-(FsPOC⁺ + FhPOC⁺), -(FsCaCO₃⁺ + FhCaCO₃⁺),
                                 -(FsbSi⁺ + FhbSi⁺), -(Fsdust⁺ + Fhdust⁺), bgc.nutrients.iron, ballast)

        @inbounds detritus.sinking_mass[i, j, k] = 𝓜

        FsCaCO₃ = FsCaCO₃⁺ * dCaCO₃ - ΠCaCO₃ * (one(γ_calcite) - γ_calcite) * (one(dCaCO₃) - dCaCO₃) * ℓCaCO₃
        FhCaCO₃ = FhCaCO₃⁺ * dhard  - ΠCaCO₃ * γ_calcite * Δz

        FsbSi = FsbSi⁺ * dbSi  - ΠbSi * (one(γ_opal) - γ_opal) * (one(dbSi) - dbSi) * ℓbSi
        FhbSi = FhbSi⁺ * dhard - ΠbSi * γ_opal * Δz

        Fsdust = Fsdust⁺ * ddust 
        Fhdust = Fhdust⁺ * dhdust

        A = ΠPOC - ρ_calcite * ΠCaCO₃ - ρ_opal * ΠbSi

        dustin  = -(Fsdust⁺ + Fhdust⁺)
        dustout = -(Fsdust + Fhdust)

        𝒬 = ifelse((𝒬 > ∅) & (dustin > ∅), 𝒬 * dustout / dustin, zero(𝒬))

        paying  = (𝒬 > ∅) & (ΠPOC > ∅)
        𝒬paid   = 𝒬 - A * Δz
        cleared = 𝒬paid < ∅

        A = ifelse(paying, ifelse(cleared, -𝒬paid * Δz⁻, zero(A)), A)
        𝒬 = ifelse(paying, max(𝒬paid, zero(𝒬paid)), 𝒬)

        FsPOC = FsPOC⁺ * dPOC - A * (one(dPOC) - dPOC) * ℓPOC

        protected = ρ_calcite * -(FsCaCO₃ + FhCaCO₃) + ρ_opal * -(FsbSi + FhbSi) + ρ_dust * dustout - 𝒬
        FhPOC     = -max(protected, zero(protected))

        FhPOC = ifelse((FhPOC⁺ == ∅) & (ΠPOC == ∅), zero(FhPOC), FhPOC)

        RCaCO₃ = ΠCaCO₃ - ((FsCaCO₃⁺ - FsCaCO₃) + (FhCaCO₃⁺ - FhCaCO₃)) * Δz⁻
        RbSi   = ΠbSi   - ((FsbSi⁺   - FsbSi)   + (FhbSi⁺   - FhbSi))   * Δz⁻
        RPOC   = ΠPOC   - ((FsPOC⁺   - FsPOC)   + (FhPOC⁺   - FhPOC))   * Δz⁻
        Rdust  =        - ((Fsdust⁺  - Fsdust)  + (Fhdust⁺  - Fhdust))  * Δz⁻

        ΠPFe = particulate_iron_production(i, j, k, grid, bgc.plankton, bgc, model_fields, aux)

        FPOCin = -(FsPOC⁺ + FhPOC⁺)
        FPFein = -FPFe⁺

        RPFe = ifelse(FPOCin > ∅,
                      RPOC * FPFein / FPOCin,
                      RPOC * ballast.iron_to_carbon_reference)

        RPFe += FPFein * ballast.iron_desorption_rate

        FPFe = FPFe⁺ - Δz * (ΠPFe - RPFe)

        clipFe = FPFe > ∅
        RPFe   = ifelse(clipFe, FPFein * Δz⁻ + ΠPFe, RPFe)
        FPFe   = ifelse(clipFe, zero(FPFe), FPFe)

        RPFe += Rdust * ballast.dust_to_iron +
                get_sedimentary_iron_flux(i, j, k, grid, ballast.sedimentary_iron_flux) * Δz⁻

        FPOPin = -FPOP⁺

        RPOP = ifelse(FPOCin > ∅,
                      RPOC * FPOPin / FPOCin,
                      ifelse(ΠPOC > ∅, RPOC * ΠPOP / ΠPOC, zero(RPOC)))

        FPOP = FPOP⁺ - Δz * (ΠPOP - RPOP)

        clipP = FPOP > ∅
        RPOP  = ifelse(clipP, FPOPin * Δz⁻ + ΠPOP, RPOP)
        FPOP  = ifelse(clipP, zero(FPOP), FPOP)

        at_floor = (k == k_bottom) & !ballast.open_bottom

        RPOC   += -(FsPOC   + FhPOC)   * Δz⁻ * at_floor
        RPOP   += -FPOP                * Δz⁻ * at_floor
        RbSi   += -(FsbSi   + FhbSi)   * Δz⁻ * at_floor
        RCaCO₃ += -(FsCaCO₃ + FhCaCO₃) * Δz⁻ * at_floor
        RPFe   += -FPFe                * Δz⁻ * at_floor

        FsPOC   *= !at_floor; FhPOC   *= !at_floor
        FsbSi   *= !at_floor; FhbSi   *= !at_floor
        FsCaCO₃ *= !at_floor; FhCaCO₃ *= !at_floor
        FPOP    *= !at_floor; FPFe    *= !at_floor

        @inbounds begin
            detritus.remineralisation.POC[i, j, k]     = RPOC
            detritus.remineralisation.POP[i, j, k]     = RPOP
            detritus.remineralisation.bSi[i, j, k]     = RbSi
            detritus.remineralisation.PFe[i, j, k]     = RPFe
            inorganic_carbon.remineralisation[i, j, k] = RCaCO₃
        end

        # no-op if not saving fluxes
        store_fluxes!(detritus.fluxes, i, j, k, FsPOC, FhPOC, FsCaCO₃ + FhCaCO₃, FsbSi + FhbSi,
                      Fsdust + Fhdust, FPFe, FPOP)
    end

    @inbounds begin
        detritus.floor_flux.POC[i, j, 1]     = -(FsPOC + FhPOC)
        detritus.floor_flux.POP[i, j, 1]     = -FPOP
        detritus.floor_flux.bSi[i, j, 1]     = -(FsbSi + FhbSi)
        detritus.floor_flux.PFe[i, j, 1]     = -FPFe
        inorganic_carbon.floor_flux[i, j, 1] = -(FsCaCO₃ + FhCaCO₃)
    end
end

function update_biogeochemical_state!(model, ::MultiElementRefractoryDissolved,
                                      npd::NutrientsPlanktonDetritus)
    grid = model.grid
    
    launch!(architecture(grid), grid, :xy, compute_ballast_fluxes!,
            grid, npd, fields(model), biogeochemical_auxiliary_fields(model.biogeochemistry))

    return nothing
end
