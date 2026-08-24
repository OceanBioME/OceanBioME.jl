#####
##### Calcite as a ballast mineral rather than a tracer
#####
##### A sibling of `ExplicitCalciumCarbonate` with the calcite tracer removed: calcite is one of the three
##### ballast minerals, so its sinking flux is solved by the ballast sweep alongside organic carbon, opal
##### and dust, and only its remineralisation reaches dissolved inorganic carbon and alkalinity.
#####
##### The saturation state is still computed every step — not for dissolution, since the length scale
##### replaces that, but because a sediment's calcite burial threshold needs it.
#####
##### ONE remineralisation field serves every realisation. The calcite flux and remineralisation are
##### provably identical across realisations, because the production, hard fraction, ballast ratio and
##### dissolution length are all shared; only the carbon and alkalinity differ.
#####

struct ImplicitExplicitCalcite{N, CC, RM, FL, SS} <: AbstractInorganicCarbon
      carbon_chemistry :: CC
      remineralisation :: RM   # filled by the ballast sweep; shared by every realisation
            floor_flux :: FL   # what reaches the sea floor, positive-down
    calcite_saturation :: SS   # per realisation; a sediment's burial threshold needs it

    ImplicitExplicitCalcite{N}(carbon_chemistry::CC, remineralisation::RM,
                               floor_flux::FL, calcite_saturation::SS) where {N, CC, RM, FL, SS} =
        new{N, CC, RM, FL, SS}(carbon_chemistry, remineralisation, floor_flux, calcite_saturation)
end

"""
    ImplicitExplicitCalcite(grid; replicates = 1, carbon_chemistry = CarbonChemistry())

The inorganic carbon component for implicit ballast sinking
([`DissolvedRefractoryImplicitCNP`](@ref)), tracking `DIC` and `Alk` only.

Unlike [`ExplicitCalciumCarbonate`](@ref) there is **no `CaCO₃` tracer**: calcite is a ballast mineral,
so its sinking flux is solved implicitly by the ballast sweep and only the resulting remineralisation
returns to `DIC`/`Alk`.

Calcifiers form calcite as they grow, removing `DIC` and twice that from `Alk`, routed through the same
three per-plankton hooks [`ExplicitCalciumCarbonate`](@ref) uses. Calcite leaving the living pool as
solid waste enters the sinking flux; calcite dissolved in zooplankton guts returns straight to
`DIC`/`Alk`.

`replicates > 1` manifests independent `DIC`/`Alk` realisations (`DIC1`, `Alk1`, `DIC2`, …).
"""
function ImplicitExplicitCalcite(grid::AbstractGrid{FT};
                        replicates = 1,
                        carbon_chemistry = CarbonChemistry(FT)) where FT

    field_names = saturation_field_names(Val(replicates))
    calcite_saturation = NamedTuple{field_names}(ntuple(_ -> CenterField(grid), replicates))

    remineralisation = CenterField(grid)
    floor_flux       = Field{Center, Center, Nothing}(grid)

    manifest_ballast_calcite_replicates!(replicates)

    return ImplicitExplicitCalcite{replicates}(carbon_chemistry, remineralisation, floor_flux, calcite_saturation)
end

const NPD_BC = NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:Any, <:ImplicitExplicitCalcite}

required_biogeochemical_tracers(::ImplicitExplicitCalcite{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::ImplicitExplicitCalcite{N}) where N =
    (map(n -> Symbol(:DIC, n), 1:N)..., map(n -> Symbol(:Alk, n), 1:N)...)

saturation_field_names(::ImplicitExplicitCalcite{N}) where N = saturation_field_names(Val(N))

carbonate_field_names(::ImplicitExplicitCalcite{1}) = ((:DIC, :Alk), )
carbonate_field_names(::ImplicitExplicitCalcite{N}) where N = ntuple(n -> (Symbol(:DIC, n), Symbol(:Alk, n)), N)

required_biogeochemical_auxiliary_fields(ic::ImplicitExplicitCalcite) =
    (saturation_field_names(ic)..., :CaCO₃_remin, :CaCO₃_floor)

biogeochemical_auxiliary_fields(ic::ImplicitExplicitCalcite) =
    merge(ic.calcite_saturation, (CaCO₃_remin = ic.remineralisation, CaCO₃_floor = ic.floor_flux))

@inline calcite_remineralisation(i, j, k, grid, ic::ImplicitExplicitCalcite, fields) =
    @inbounds ic.remineralisation[i, j, k]

@inline net_calcium_carbonate_production(i, j, k, grid, bgc::NPD_BC, fields, auxiliary_fields) =
    biological_calcium_carbonate_precipitation(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) -
    calcite_remineralisation(i, j, k, grid, bgc.inorganic_carbon, fields) -
    biological_calcium_carbonate_dissolution(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

function update_biogeochemical_state!(model, ic::ImplicitExplicitCalcite, npd::NutrientsPlanktonDetritus)
    for (saturation, carbonate_names) in zip(values(ic.calcite_saturation), carbonate_field_names(ic))
        compute_calcium_carbonate_saturation!(ic.carbon_chemistry, saturation, model, carbonate_names)
    end

    return nothing
end

const _manifested_ballast_calcite = Set{Int}()

function manifest_ballast_calcite_replicates!(N)
    (N > 1 && !(N in _manifested_ballast_calcite)) || return nothing
    push!(_manifested_ballast_calcite, N)

    for n in 1:N
        DIC_name = Symbol(:DIC, n)
        Alk_name = Symbol(:Alk, n)

        @eval begin
            @inline (bgc::NPD_BC)(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) -
                net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields)

            @inline (bgc::NPD_BC)(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields) -
                2 * net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields)
        end
    end

    return nothing
end

Adapt.adapt_structure(to, ic::ImplicitExplicitCalcite{N}) where N =
    ImplicitExplicitCalcite{N}(adapt(to, ic.carbon_chemistry),
                      adapt(to, ic.remineralisation),
                      adapt(to, ic.floor_flux),
                      adapt(to, ic.calcite_saturation))

Base.summary(ic::ImplicitExplicitCalcite{1}) = string("ImplicitExplicitCalcite $(required_biogeochemical_tracers(ic))")
Base.summary(ic::ImplicitExplicitCalcite{N}) where N =
    string("ImplicitExplicitCalcite{realisations = $N} $(required_biogeochemical_tracers(ic))")

function Base.show(io::IO, ic::ImplicitExplicitCalcite{N}) where N
    msg = "ImplicitExplicitCalcite $(required_biogeochemical_tracers(ic)) (calcite implicit — a ballast flux)"

    N > 1 && (msg *= "\n└── Realisations: $N")

    print(io, msg)

    return nothing
end

group_element_tracers(::ImplicitExplicitCalcite, args...) = NamedTuple()

group_element_tracers(::ImplicitExplicitCalcite{1}, bgc::NPD{FT}, ::Val{:carbon}) where FT = (; DIC = one(FT))

function group_element_tracers(ic::ImplicitExplicitCalcite{N}, bgc::NPD{FT}, ::Val{:carbon}) where {N, FT}
    names = filter(n -> startswith(String(n), "DIC"), required_biogeochemical_tracers(ic))

    return NamedTuple{names}(ntuple(_ -> one(FT), length(names)))
end
