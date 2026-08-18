using OceanBioME, Oceananigans, Test

using OceanBioME: conserved_tracers
using OceanBioME.Models.NutrientsPlanktonDetritusModels: InstantRemineralisationDetritus

const NPDModels = OceanBioME.Models.NutrientsPlanktonDetritusModels
const DetritusModels = NPDModels.DetritusModels

struct ExternalContractPlankton end

Oceananigans.Biogeochemistry.required_biogeochemical_tracers(::ExternalContractPlankton) =
    (:P1, :P2, :C1, :C2)

Oceananigans.Biogeochemistry.required_biogeochemical_auxiliary_fields(::ExternalContractPlankton) = tuple()

NPDModels.group_element_tracers(::ExternalContractPlankton, bgc, ::Val{:nitrogen}) =
    (P1 = 1.0, P2 = 1.0, C1 = 2.0, C2 = 2.0)

NPDModels.chlorophyll_ratio(::ExternalContractPlankton) = 1.31
NPDModels.chlorophyll_ratio(::ExternalContractPlankton, ::Val{:P2}) = 0.9

const ExternalContractNPD{FT} = NPDModels.NutrientsPlanktonDetritus{FT, <:Any, <:ExternalContractPlankton}
const ExternalContractTracer = Union{Val{:P1}, Val{:P2}, Val{:C1}, Val{:C2}}

@inline (bgc::ExternalContractNPD{FT})(i, j, k, grid, ::ExternalContractTracer, clock, fields, auxiliary_fields) where FT = FT(1)

@inline NPDModels.nutrient_uptake(i, j, k, grid, plankton::ExternalContractPlankton, bgc::ExternalContractNPD{FT}, fields, auxiliary_fields) where FT = FT(0.4)
@inline NPDModels.solid_waste(i, j, k, grid, plankton::ExternalContractPlankton, bgc::ExternalContractNPD{FT}, fields, auxiliary_fields) where FT = FT(0.3)
@inline NPDModels.dissolved_waste(i, j, k, grid, plankton::ExternalContractPlankton, bgc::ExternalContractNPD{FT}, fields, auxiliary_fields) where FT = FT(0.2)
@inline NPDModels.inorganic_waste(i, j, k, grid, plankton::ExternalContractPlankton, bgc::ExternalContractNPD{FT}, fields, auxiliary_fields) where FT = FT(0.1)

@inline DetritusModels.grazing(i, j, k, grid, ::Val{:sPOM}, plankton::ExternalContractPlankton, bgc::ExternalContractNPD{FT}, fields, auxiliary_fields) where FT = FT(0.5)

@testset "External plankton interface" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    plankton = ExternalContractPlankton()
    bgc = NPDModels.NutrientsPlanktonDetritus{Float64}(
        Nutrients(; nitrogen = OceanBioME.Models.N),
        plankton,
        InstantRemineralisationDetritus(),
        nothing,
        nothing,
    )

    fields = NamedTuple()
    auxiliary_fields = NamedTuple()

    @test Oceananigans.Biogeochemistry.required_biogeochemical_tracers(bgc) == (:N, :P1, :P2, :C1, :C2)
    @test bgc(1, 1, 1, grid, Val(:P1), nothing, fields, auxiliary_fields) == 1.0
    @test NPDModels.nutrient_uptake(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields) == 0.4
    @test NPDModels.solid_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields) == 0.3
    @test NPDModels.dissolved_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields) == 0.2
    @test NPDModels.inorganic_waste(1, 1, 1, grid, plankton, bgc, fields, auxiliary_fields) == 0.1
    @test DetritusModels.grazing(1, 1, 1, grid, Val(:sPOM), plankton, bgc, fields, auxiliary_fields) == 0.5
    @test NPDModels.chlorophyll_ratio(plankton, Val(:P1)) == 1.31
    @test NPDModels.chlorophyll_ratio(plankton, Val(:P2)) == 0.9
end

@testset "External plankton conservation metadata" begin
    nutrients = Nutrients(; nitrogen = OceanBioME.Models.N,
                           phosphate = OceanBioME.Models.PO₄,
                           iron = OceanBioME.Models.Fe,
                           silicate = OceanBioME.Models.Si)

    bgc = NPDModels.NutrientsPlanktonDetritus{Float64}(
        nutrients,
        ExternalContractPlankton(),
        InstantRemineralisationDetritus(),
        CarbonateSystem(),
        Oxygen(),
    )

    tracer_groups = conserved_tracers(bgc)

    @test (tracer_groups.nitrogen.P1, tracer_groups.nitrogen.P2,
           tracer_groups.nitrogen.C1, tracer_groups.nitrogen.C2) == (1.0, 1.0, 2.0, 2.0)
    @test (tracer_groups.phosphate.P1, tracer_groups.iron.P1,
           tracer_groups.silicate.P1, tracer_groups.carbon.P1) == (1 / 16, 0.0032 / 16, 0.0, 106 / 16)
end
