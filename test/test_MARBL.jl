include("dependencies_for_runtests.jl")

using OceanBioME.Models.NutrientsPlanktonDetritusModels: required_biogeochemical_tracers
using Adapt: adapt
using Oceananigans.Fields: CenterField

const DET = OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels.MultiElement

# a rate which is non-linear in light, and one which is linear, for the sub column tests below
@subcolumn_average light_saturation(PAR) = 1 - exp(-PAR)
@subcolumn_average scaled_light(PAR, a) = a * PAR

#####
##### A closed box, so every element must be conserved exactly. Sinking is off and the bottom is closed.
#####

const BOX_CC = CarbonChemistry(; density_function = (args...) -> 1026.0)

function marbl_box(; PAR = 100.0, plankton = nothing, grid = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1)))
    plankton = isnothing(plankton) ? ManyPhytoZoo() : plankton

    # explicitly the EXPLICIT-sinking configuration: the default `MARBL` sweeps a column for the
    # implicit ballast flux, which a one cell box is no home for (and whose particulates are not tracers)
    bgc = MARBL_ExplicitSinking(grid;
                plankton,
                sinking_speed = 0.0,
                nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, Fe, Si),
                inorganic_carbon = ExplicitCalciumCarbonate(grid; carbon_chemistry = BOX_CC),
                light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(PAR)))

    return NonhydrostaticModel(grid; biogeochemistry = bgc, advection = nothing, buoyancy = nothing,
                               tracers = (:T, ),
                               auxiliary_fields = (S = ConstantField(35.0), ))
end

# the initial state used by every conservation test
const BOX_STATE = (T = 12.0, NO₃ = 5.0, NH₄ = 0.2, PO₄ = 0.5, Fe = 2e-4, Si = 3.0,
                   spC = 0.6, spChl = 0.15, spP = 0.005, spFe = 2e-5, spCaCO₃ = 0.02,
                   diatC = 0.4, diatChl = 0.12, diatP = 0.004, diatFe = 1.5e-5, diatSi = 0.05,
                   diazC = 0.02, diazChl = 0.004, diazP = 0.0002, diazFe = 1e-6, zooC = 0.3,
                   DOC = 30.0, DON = 4.0, DOP = 0.3, DOCr = 40.0, DONr = 2.0, DOPr = 0.1,
                   POC = 0.5, POP = 0.004, bSi = 0.2, PFe = 1e-5,
                   DIC = 2000.0, Alk = 2300.0, CaCO₃ = 0.02, O₂ = 250.0)

set_box!(model, state = BOX_STATE) =
    for name in keys(model.tracers)
        haskey(state, name) && set!(model.tracers[name], state[name])
    end

# element sums, written out explicitly rather than through `conserved_tracers`, because the quotas are
# per-PFT and variable so a fixed-coefficient sum cannot express them
value(model, name) = @inbounds model.tracers[name][1, 1, 1]

carbon(m) = sum(value(m, n) for n in (:spC, :diatC, :diazC, :zooC, :DOC, :DOCr, :POC, :DIC)) +
            value(m, :spCaCO₃) + value(m, :CaCO₃)

phosphorus(m) = sum(value(m, n) for n in (:spP, :diatP, :diazP, :DOP, :DOPr, :POP, :PO₄)) +
                value(m, :zooC) / 117

iron(m) = sum(value(m, n) for n in (:spFe, :diatFe, :diazFe, :PFe, :Fe)) +
          3.0e-6 * value(m, :zooC)

silicon(m) = value(m, :diatSi) + value(m, :bSi) + value(m, :Si)

@testset "MARBL" begin
    @testset "construction and tracer lists" begin
        grid = RectilinearGrid(size = (1, 1, 4), extent = (1, 1, 100))

        bgc = MARBL(grid)
        tracers = required_biogeochemical_tracers(bgc)

        # the plankton tracers follow each PFT's own traits
        for name in (:spC, :spChl, :spP, :spFe, :spCaCO₃, :diatSi, :diazC, :zooC)
            @test name in tracers
        end

        @test !(:diazSi in tracers)      # only the silicifier carries silicon
        @test !(:diatCaCO₃ in tracers)   # only the calcifier carries calcite
        @test !(:T in tracers)           # temperature comes from the physics

        # the iron cycle and the redox chemistry are on by default
        @test :Lig in tracers
        @test :O₂ in tracers

        # sinking is MARBL's implicit ballast flux by default, so the particulates and the detrital
        # calcite are diagnostic, not tracers
        for name in (:POC, :POP, :bSi, :PFe, :CaCO₃)
            @test !(name in tracers)
        end

        # ... while the explicitly sinking testing configuration carries them
        explicit_tracers = required_biogeochemical_tracers(MARBL_ExplicitSinking(grid))

        for name in (:POC, :POP, :bSi, :PFe, :CaCO₃)
            @test name in explicit_tracers
        end

        # the +cocco variant: coccolithophores calcify and the small phytoplankton do not
        cocco = MARBL_Cocco(grid)
        cocco_tracers = required_biogeochemical_tracers(cocco)

        @test :coccoCaCO₃ in cocco_tracers
        @test !(:spCaCO₃ in cocco_tracers)
        @test :coccoC in cocco_tracers
    end

    @testset "the model builds and steps" begin
        model = marbl_box()
        set_box!(model)

        @test model isa NonhydrostaticModel

        time_step!(model, 100.0)

        @test all(isfinite, interior(model.tracers.spC))
        @test all(isfinite, interior(model.tracers.DIC))
        @test all(isfinite, interior(model.tracers.O₂))
    end

    @testset "conservation in a closed box" begin
        model = marbl_box()
        set_box!(model)

        C₀, P₀, Fe₀, Si₀ = carbon(model), phosphorus(model), iron(model), silicon(model)

        for _ in 1:20
            time_step!(model, 100.0)
        end

        @test carbon(model)     ≈ C₀  rtol = 1e-10
        @test phosphorus(model) ≈ P₀  rtol = 1e-10
        @test iron(model)       ≈ Fe₀ rtol = 1e-10
        @test silicon(model)    ≈ Si₀ rtol = 1e-10
    end

    @testset "rate laws" begin
        MP = OceanBioME.Models.NutrientsPlanktonDetritusModels.PlanktonModels.MARBLPlankton

        # temperature scaling is a plain Q10 at the reference temperature
        @test MP.temperature_scaling(1.7, 30.0, 30.0) ≈ 1
        @test MP.temperature_scaling(1.7, 30.0, 40.0) ≈ 1.7

        # the light-limited rate saturates at the maximum and vanishes in the dark
        @test MP.light_limited_rate(1.0, 0.39, 0.02, 0.0, 0.0) ≈ 0
        @test MP.light_limited_rate(1.0, 0.39, 0.02, 1e12, 0.0) ≈ 1

        # Liebig limitation takes the minimum, and a nitrogen fixer is never nitrogen limited
        f = MP.minimum_limitation(false, false, 0.25, 0.01, 0.01, 0.3, 3e-5, 0.0,
                                  5.0, 0.2, 0.5, 0.1, 2e-4, 3.0)
        @test 0 < f < 1

        fixer = MP.minimum_limitation(false, true, 0.25, 0.01, 0.01, 0.3, 3e-5, 0.0,
                                      0.0, 0.0, 0.5, 0.1, 2e-4, 3.0)
        @test fixer > 0   # no nitrogen at all, but the fixer still grows

        # the growth iron quota economises below the threshold and is flat above it
        @test MP.growth_iron_quota(30e-6, 2.5e-6, 10.0, 3e-5, 1.0) ≈ 30e-6
        @test MP.growth_iron_quota(30e-6, 2.5e-6, 10.0, 3e-5, 0.0) ≈ 2.5e-6

        # light which is not split runs the rate once on the plain value; split light runs it in each
        # sub column and area weights the result
        grid = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1))
        field = CenterField(grid)
        set!(field, 3.0)

        plain = (PAR = field, )
        split = (PAR = SubcolumnPAR(field, (0.5, 0.5), (2.0, 0.0)), )

        saturating(aux) = light_saturation(@preserve_subcolumns(aux.PAR[1, 1, 1]))
        linear(aux)     = scaled_light(@preserve_subcolumns(aux.PAR[1, 1, 1]), 2.0)

        @test saturating(plain) ≈ 1 - exp(-3.0)
        @test saturating(split) ≈ 0.5 * (1 - exp(-6.0)) + 0.5 * (1 - exp(0.0))

        # Σⱼ φⱼrⱼ = 1, so splitting leaves a rate which IS linear in light unchanged
        @test linear(plain) ≈ 6
        @test linear(split) ≈ 6

        # and the split light still indexes like the field it wraps, giving the mean
        @test split.PAR[1, 1, 1] == 3.0
    end

    @testset "redox rate laws" begin
        OX = OceanBioME.Models.NutrientsPlanktonDetritusModels.OxygenModels

        # nitrification is off in bright water and full in the dark
        @test OX.nitrification_taper(10.0, 5.0, 1.0) ≈ 0     # both interfaces above the limit
        @test OX.nitrification_taper(0.5, 0.1, 1.0) ≈ 1      # both below
        @test 0 < OX.nitrification_taper(2.0, 0.5, 1.0) < 1  # straddling

        # the oxygen fractions are complementary and clamped
        @test OX.suboxic_fraction(100.0, 5.0, 5.0) ≈ 0
        @test OX.suboxic_fraction(0.0, 5.0, 5.0)   ≈ 1
        @test OX.aerobic_fraction(100.0, 5.0, 5.0) ≈ 1
        @test OX.aerobic_fraction(0.0, 5.0, 5.0)   ≈ 0
    end

    @testset "iron speciation" begin
        NU = OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels

        ci = ComplexedIron()
        K  = ci.ligand_stability_constant

        # no iron means no iron of either form
        @test all(NU.iron_speciation(0.0, 1e-3, K) .≈ 0)

        # speciation conserves total iron
        for (Fe, Lig) in ((2e-4, 1e-3), (5e-3, 1e-3), (1e-6, 6e-4))
            Fefree, FeLig = NU.iron_speciation(Fe, Lig, K)

            @test Fefree ≥ 0
            @test FeLig ≥ 0
            @test Fefree + FeLig ≈ Fe rtol = 1e-8
            @test FeLig ≤ Lig * (1 + 1e-8)   # bound iron cannot exceed the ligand
        end
    end

    @testset "sediment burial" begin
        s = BurialDenitrification()
        SD = OceanBioME.Models.SedimentModels

        # burial efficiency rises with the flux and is capped
        f_low  = SD.organic_carbon_burial_fraction(s, 1e-8)
        f_high = SD.organic_carbon_burial_fraction(s, 1e-3)

        @test f_low < f_high
        @test f_high ≤ s.maximum_organic_burial_fraction

        # calcite is buried only above the saturation threshold
        @test SD.calcite_burial_fraction(s, 1.5)
        @test !SD.calcite_burial_fraction(s, 0.5)

        # all the particulate iron is removed
        @test SD.iron_burial_fraction(s, 1e-9) == 1

        @test :buried_C in SD.required_sediment_fields(s)
        @test :POC in SD.sinking_fluxes(s)
    end

    @testset "ballast sinking" begin
        DM = OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels.MultiElement

        b = DET.Ballast()

        # `scale_length` takes a z-coordinate (negative below the surface); the depth factor
        # is flat outside the given range and monotonic within it
        @test DM.scale_length(0.0, b)     ≈ b.scale_length_values[1]
        @test DM.scale_length(-5000.0, b) ≈ b.scale_length_values[end]
        @test DM.scale_length(-400.0, b) > DM.scale_length(-150.0, b)

        # the low-oxygen stretch lengthens the scale and is inert in oxygenated water
        @test DM.oxygen_scale_factor(300.0, b) ≈ 1
        @test DM.oxygen_scale_factor(0.0, b)   ≈ b.oxygen_scaling_factor
        @test DM.oxygen_scale_factor(25.0, b)  > 1
    end

    @testset "adapt, isbits and show" begin
        grid = RectilinearGrid(size = (1, 1, 4), extent = (1, 1, 100))

        for component in (ManyPhytoZoo(), MARBL_small_phyto(), MARBL_zooplankton().zoo,
                          ComplexedIron(), RedoxOxygen(), BurialDenitrification(), DET.Ballast())
            @test sprint(show, component) isa String
            @test adapt(Array, component) isa Any
        end

        @test isbits(MARBL_small_phyto())
        @test isbits(ComplexedIron())
        @test isbits(RedoxOxygen())
        @test isbits(BurialDenitrification())

        @test sprint(show, MultiElementRefractoryDissolvedParticulate(grid)) isa String
    end
end
