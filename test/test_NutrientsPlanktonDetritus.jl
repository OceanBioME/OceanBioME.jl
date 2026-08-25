include("dependencies_for_runtests.jl")

using OceanBioME: conserved_tracers
using OceanBioME.Models.NutrientsPlanktonDetritusModels: InstantRemineralisationDetritus, 
                                                         CarbonNitrogenDissolvedParticulate
using Oceananigans, CUDA, Random

Random.seed!(42)

grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

test_models = []

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

using CUDA
using OceanBioME: conserved_tracers

function get_conservation_values(bgc, fields)
    tracer_groups = conserved_tracers(bgc)

    conserved_values = NamedTuple()

    for (group_name, group) in pairs(tracer_groups)
        value = 0.0

        for (name, sf) in pairs(group)
            value += sf * (CUDA.@allowscalar first(Field(Integral(fields[name]))))
        end

        conserved_values = merge(conserved_values, NamedTuple{(group_name,)}((value, )))
    end

    return conserved_values
end

# everything should be fine as long as everythings positive right???
set_default!(model) = [set!(tracer, rand()*10) for tracer in model.tracers]

function check_stepped(model)
    set_default!(model)
    initial_values = NamedTuple(name => CUDA.@allowscalar model.tracers[name][1, 1, 1]
                                for name in keys(model.tracers))
    time_step!(model, 1)
    for (name, initial) in pairs(initial_values)
        if !iszero(initial)
            @test CUDA.@allowscalar model.tracers[name][1, 1, 1] != initial
        end
    end

    return nothing
end

function check_conservations(model, n_timesteps = 100)
    set_default!(model)

    initial_values = get_conservation_values(model.biogeochemistry, model.tracers)

    for _ in 1:n_timesteps
        time_step!(model, 1)
    end

    final_values = get_conservation_values(model.biogeochemistry, model.tracers)

    for (name, value) in pairs(initial_values)
        @info "Checking conservation of $(name)"
        @test isapprox(final_values[name], value, atol = n_timesteps*eps(value))
    end
end

# need to do at least one test of detritus and dissolved particulate with advection on
nutrients_options = (Nutrients(; nitrogen = OceanBioME.N),
                     Nutrients(; nitrogen = NitrateAmmonia()),
                     Nutrients(; phosphate = OceanBioME.PO₄),
                     Nutrients(; iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = OceanBioME.N, phosphate = OceanBioME.PO₄),
                     Nutrients(; nitrogen = OceanBioME.N, iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = NitrateAmmonia(), phosphate = OceanBioME.PO₄),
                     Nutrients(; nitrogen = NitrateAmmonia(), iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = OceanBioME.N, phosphate = OceanBioME.PO₄, iron = OceanBioME.Fe),
                     Nutrients(; nitrogen = NitrateAmmonia(), phosphate = OceanBioME.PO₄, iron = OceanBioME.Fe))

detritus_options = (InstantRemineralisationDetritus(),
                    Detritus(grid), 
                    DissolvedParticulate(grid, :DOP, :POP), 
                    DissolvedParticulate(grid),
                    CarbonNitrogenDissolvedParticulate(grid))

# TODO: test multiple separately
inorganic_carbon_options = (nothing, 
                            CarbonateSystem(),)#, CarbonateSystem(2))

oxygen_options = (nothing, 
                  Oxygen(),)

@testset "Elemental conservations" begin 
    # maybe this is pointless but would like to keep Abiotic working for testing
    # (abiotic is going to result in everything being remineralised since theres no path to make detritus)
    for nutrients in nutrients_options[end-3:end], # empty slot, nitrate/ammonia, and N
        detritus in detritus_options,
        inorganic_carbon in inorganic_carbon_options,
        oxygen in oxygen_options

        biogeochemistry = NutrientsPlanktonDetritus(grid;
                                                    plankton = Abiotic(),
                                                    nutrients,
                                                    detritus,
                                                    inorganic_carbon,
                                                    oxygen)

        @info summary(biogeochemistry.underlying_biogeochemistry)

        model = NonhydrostaticModel(grid; 
                                    advection = nothing,
                                    biogeochemistry)

        time_step!(model, 1.0) 
        @test (CUDA.@allowscalar all([all(values .== 0) for values in values(model.tracers)]))

        check_conservations(model)
    end

    for plankton in (PhytoZoo, ImplicitProductivity),
        nutrients in nutrients_options[end-3:end], # empty slot, nitrate/ammonia, and N
        detritus in detritus_options,
        inorganic_carbon in inorganic_carbon_options,
        oxygen in oxygen_options

        nutrient_half_saturations = 
            NamedTuple(name=>half_sat for 
                        (name, half_sat) in pairs((phosphate = 1, iron = 0.001)) 
                        if !isnothing(getproperty(nutrients, name)))

        if nutrients.nitrogen isa NitrateAmmonia
            nutrient_half_saturations = merge(nutrient_half_saturations, (nitrate = 0.5, ammonia = 0.1))
        elseif nutrients.nitrogen isa OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels.SingleTracerNutrient
            nutrient_half_saturations = merge(nutrient_half_saturations, (; nitrate = 0.5))
        end

        biogeochemistry = NutrientsPlanktonDetritus(grid;
                                                    plankton = plankton(; nutrient_half_saturations),
                                                    nutrients,
                                                    detritus,
                                                    inorganic_carbon,
                                                    oxygen,
                                                    light_attenuation)

        @info summary(biogeochemistry.underlying_biogeochemistry)

        model = NonhydrostaticModel(grid; 
                                    advection = nothing,
                                    biogeochemistry)

        time_step!(model, 1.0) 
        @test (CUDA.@allowscalar all([all(values .== 0) for values in values(model.tracers)]))

        check_conservations(model)

        if !((detritus isa InstantRemineralisationDetritus)&(plankton == ImplicitProductivity))
            check_stepped(model)
        end
    end
end

@testset "Alkalinity sign" begin
    biogeochemistry = LOBSTER(grid; inorganic_carbon = CarbonateSystem(), plankton = PhytoZoo(rain_ratio = 0))

    model = NonhydrostaticModel(grid; advection = nothing, biogeochemistry)

    set!(model, NO₃ = 10, NH₄ = 0, P = 1, Z = 0.1, DOM = 0.1, sPOM = 0.1, bPOM = 0.1, DIC = 2000, Alk = 1900)

    NO₃₀, NH₄₀, Alk₀ = CUDA.@allowscalar model.tracers.NO₃[1, 1, 1], model.tracers.NH₄[1, 1, 1], model.tracers.Alk[1, 1, 1]
    time_step!(model, 1) 
    ΔNO₃, ΔNH₄, ΔAlk = CUDA.@allowscalar model.tracers.NO₃[1, 1, 1] - NO₃₀, model.tracers.NH₄[1, 1, 1] - NH₄₀, model.tracers.Alk[1, 1, 1] - Alk₀

    @test ΔNO₃ < 0 && ΔAlk > 0 # net nitrate uptake raises alkalinity
    @test isapprox(ΔAlk, -(1 + 1/16) * ΔNO₃ + (1 - 1/16) * ΔNH₄, atol=1e-12)
end

@testset "Explicit particle sinking" begin # not sure how essential this is but seems worth doing
    grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 2))

    biogeochemistry = NPZD(grid; 
                           detritus = Detritus(grid;
                                               sinking_speed=1,
                                               remineralisation_rate=0.0,
                                               open_bottom = true))

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, D = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.D[1, 1, 2], exp(-0.01*100), atol = 1e-3)
    end 

    biogeochemistry = NPZD(grid; 
                           detritus = Detritus(grid;
                                               sinking_speed=1,
                                               remineralisation_rate=0.0,
                                               open_bottom = false))

    # we scale to zero at the bottom which in a normal situation is fine but in this situation means the top cell is also not as set
    first_cell_flux_velocity = CUDA.@allowscalar biogeochemistry.underlying_biogeochemistry.detritus.sinking_speeds.w[1, 1, 2]

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, D = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.D[1, 1, 2], exp(0.01*100 * first_cell_flux_velocity), atol = 1e-3)
        @test Field(Integral(model.tracers.D))[1, 1, 1] ≈ 1
    end  

    biogeochemistry = NPZD(grid; 
                           detritus = DissolvedParticulate(grid;
                                                           sinking_speeds=(0.1, 1),
                                                           particulate_remineralisation_rate=(0.0, 0.0),
                                                           open_bottom = true))

    model = NonhydrostaticModel(grid; 
                                biogeochemistry, 
                                advection = UpwindBiased(order=1),
                                timestepper = :QuasiAdamsBashforth2)

    set!(model, sPOM = (x, y, z) -> z>-1, bPOM = (x, y, z) -> z>-1)

    for _ in 1:100
        time_step!(model, 0.01)
    end

    CUDA.@allowscalar begin
        @test isapprox(model.tracers.sPOM[1, 1, 2], exp(-0.01*100*0.1), atol = 1e-3)
        @test isapprox(model.tracers.bPOM[1, 1, 2], exp(-0.01*100), atol = 1e-3)
    end 
end

using OceanBioME.Models.NutrientsPlanktonDetritusModels: SingleTracerNutrient

@testset "NutrientsPlanktonDetritus constructors and float types" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

    lobster = LOBSTER(grid)
    npzd = NPZD(grid)
    implicit = ImplicitBiology(grid)

    @test lobster isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test npzd isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test implicit isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}

    @test lobster.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:NitrateAmmonia, Nothing, Nothing, Nothing}, <:PhytoZoo, <:DissolvedParticulate{1, 2}}
    @test npzd.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:SingleTracerNutrient, Nothing, Nothing, Nothing}, <:PhytoZoo, <:Detritus}
    @test implicit.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float64, <:Nutrients{<:SingleTracerNutrient, <:SingleTracerNutrient, <:SingleTracerNutrient, Nothing}, <:ImplicitProductivity, <:DissolvedParticulate{1, 1}}

    grid = RectilinearGrid(architecture, Float32; size=(1, 1, 1), extent=(1, 1, 2))

    lobster = LOBSTER(grid; oxygen = Oxygen(Float32))
    npzd = NPZD(grid)
    implicit = ImplicitBiology(grid)

    @test lobster isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test npzd isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}
    @test implicit isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsPlanktonDetritus}

    @test lobster.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:NitrateAmmonia, Nothing, Nothing, Nothing}, <:PhytoZoo, <:DissolvedParticulate{1, 2}}
    @test npzd.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:SingleTracerNutrient, Nothing, Nothing, Nothing}, <:PhytoZoo, <:Detritus}
    @test implicit.underlying_biogeochemistry isa NutrientsPlanktonDetritus{Float32, <:Nutrients{<:SingleTracerNutrient, <:SingleTracerNutrient, <:SingleTracerNutrient, Nothing}, <:ImplicitProductivity, <:DissolvedParticulate{1, 1}}

    # light is converted
    par = lobster.light_attenuation
    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32

    # nutrients is converted
    @test lobster.underlying_biogeochemistry.nutrients.nitrogen |> ((::NitrateAmmonia{FT}) where FT) -> FT == Float32

    # plankton is converted
    @test lobster.underlying_biogeochemistry.plankton |> ((::PhytoZoo{<:Any, <:Any, <:Any, FT}) where FT) -> FT == Float32
    @test implicit.underlying_biogeochemistry.plankton |> ((::ImplicitProductivity{FT}) where FT) -> FT == Float32

    # detritus
    @test npzd.underlying_biogeochemistry.detritus |> ((::Detritus{FT}) where FT) -> FT == Float32
    @test lobster.underlying_biogeochemistry.detritus |> ((::DissolvedParticulate{<:Any, <:Any, <:Any, <:Any, FT}) where FT) -> FT == Float32

    # oxygen (carbonate system doesn't carry any numbers)
    @test lobster.underlying_biogeochemistry.oxygen |> ((::Oxygen{FT}) where FT) -> FT == Float32
end

using Oceananigans.Units: day
using Oceananigans.Grids: znode, Center
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields, update_biogeochemical_state!
using OceanBioME.Models.CarbonChemistryModel: calcium_carbonate_saturation

# a single-cell explicit-calcium-carbonate model with prescribed T/S (T and S are required for the carbon
# chemistry). Ω is not automatically updated, so tests that call the tendency directly set it by hand.
function explicit_calcium_carbonate_model(grid;
                                plankton = Abiotic(),
                                nutrients = Nutrients(nothing, nothing, nothing, nothing),
                                detritus = InstantRemineralisationDetritus(),
                                T = 15.0, S = 35.0,
                                kwargs...)
    inorganic_carbon = ExplicitCalciumCarbonate(grid; kwargs...)

    biogeochemistry = NutrientsPlanktonDetritus(grid;
        plankton, nutrients, detritus, inorganic_carbon,
        light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100)))

    return NonhydrostaticModel(grid;
        biogeochemistry, advection = nothing,
        auxiliary_fields = (T = ConstantField(T), S = ConstantField(S)))
end

@testset "ExplicitCalciumCarbonate" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

    @testset "construction and tracers" begin
        model  = explicit_calcium_carbonate_model(grid)
        model2 = explicit_calcium_carbonate_model(grid; replicates = 2)

        @test :CaCO₃ in keys(model.tracers)
        @test (:DIC, :Alk, :CaCO₃) ⊆ keys(model.tracers)

        @test (:DIC1, :DIC2, :Alk1, :Alk2, :CaCO₃1, :CaCO₃2) ⊆ keys(model2.tracers)
        @test model2.biogeochemistry.underlying_biogeochemistry.inorganic_carbon isa ExplicitCalciumCarbonate{2}
    end

    @testset "dissolution rate law (machine precision)" begin
        k, m = 0.2 / day, 1.5
        model = explicit_calcium_carbonate_model(grid;
            calcium_carbonate_dissolution_rate = k, calcium_carbonate_dissolution_exponent = m)

        set!(model, DIC = 2000, Alk = 2300, CaCO₃ = 3.0)

        bgc = model.biogeochemistry.underlying_biogeochemistry
        Ωfield = bgc.inorganic_carbon.calcium_carbonate_saturation[1]

        Ωval  = 0.4                # undersaturated → dissolution active
        CaCO₃ = 3.0
        set!(Ωfield, Ωval)

        aux = biogeochemical_auxiliary_fields(model.biogeochemistry)
        expected_dissolution = k * max(0, 1 - Ωval)^m * CaCO₃

        tCaCO₃ = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃), model.clock, model.tracers, aux)
        tDIC   = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:DIC),   model.clock, model.tracers, aux)
        tAlk   = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:Alk),   model.clock, model.tracers, aux)

        # calcium carbonate dissolves; DIC and Alk gain 1× and 2× that (no biology ⇒ no other carbon/alkalinity source)
        @test tCaCO₃ ≈ -expected_dissolution
        @test tDIC   ≈  expected_dissolution
        @test tAlk   ≈  2 * expected_dissolution

        # supersaturated ⇒ no dissolution and, by default, no precipitation
        set!(Ωfield, 2.5)
        @test (CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃), model.clock, model.tracers, aux)) == 0
    end

    @testset "abiotic precipitation rate law" begin
        k, n = 0.1 / day, 2.0
        model = explicit_calcium_carbonate_model(grid;
            calcium_carbonate_precipitation_rate = k, calcium_carbonate_precipitation_exponent = n)

        set!(model, DIC = 2000, Alk = 2300, CaCO₃ = 3.0)

        bgc = model.biogeochemistry.underlying_biogeochemistry
        Ωfield = bgc.inorganic_carbon.calcium_carbonate_saturation[1]

        Ωval  = 2.5                # supersaturated → precipitation active
        CaCO₃ = 3.0
        set!(Ωfield, Ωval)

        aux = biogeochemical_auxiliary_fields(model.biogeochemistry)
        expected_precipitation = k * max(0, Ωval - 1)^n * CaCO₃

        tCaCO₃ = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃), model.clock, model.tracers, aux)
        tDIC   = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:DIC),   model.clock, model.tracers, aux)
        tAlk   = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:Alk),   model.clock, model.tracers, aux)

        @test tCaCO₃ ≈  expected_precipitation
        @test tDIC   ≈ -expected_precipitation
        @test tAlk   ≈ -2 * expected_precipitation
    end

    @testset "Ω matches a direct carbon-chemistry call" begin
        model = explicit_calcium_carbonate_model(grid; T = 12.0, S = 34.0)
        set!(model, DIC = 2100, Alk = 2350, CaCO₃ = 4.0)

        update_biogeochemical_state!(model.biogeochemistry, model)

        bgc = model.biogeochemistry.underlying_biogeochemistry
        Ωfield = bgc.inorganic_carbon.calcium_carbonate_saturation[1]

        z = znode(1, 1, 1, grid, Center(), Center(), Center())
        P = abs(z) * Oceananigans.defaults.gravitational_acceleration * 1026 / 100000
        Ωdirect = calcium_carbonate_saturation(bgc.inorganic_carbon.carbon_chemistry;
                                     DIC = 2100.0, T = 12.0, S = 34.0, Alk = 2350.0, P)

        @test (CUDA.@allowscalar Ωfield[1, 1, 1]) ≈ Ωdirect
    end

    @testset "biological calcium carbonate production (PhytoZoo, machine precision)" begin
        NPDM = OceanBioME.Models.NutrientsPlanktonDetritusModels
        function pz_model(ρ)
            biogeochemistry = NutrientsPlanktonDetritus(grid;
                plankton = PhytoZoo(grid; rain_ratio = ρ),
                nutrients = Nutrients(; nitrogen = OceanBioME.N),
                detritus = DissolvedParticulate(grid),
                inorganic_carbon = ExplicitCalciumCarbonate(grid),      # precipitation off by default
                light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100)))
            m = NonhydrostaticModel(grid; biogeochemistry, advection = nothing,
                auxiliary_fields = (T = ConstantField(15.0), S = ConstantField(35.0)))
            set!(m, N = 5, P = 1, Z = 0.5, DOM = 1, sPOM = 1, bPOM = 1, DIC = 2000, Alk = 2300, CaCO₃ = 3)
            return m
        end

        ρ = 0.1
        model_ρ = pz_model(ρ)
        model_0 = pz_model(0.0)

        bgc  = model_ρ.biogeochemistry.underlying_biogeochemistry
        bgc0 = model_0.biogeochemistry.underlying_biogeochemistry
        pl   = bgc.plankton

        Ωfield_ρ = bgc.inorganic_carbon.calcium_carbonate_saturation[1]
        Ωfield_0 = bgc0.inorganic_carbon.calcium_carbonate_saturation[1]
        set!(Ωfield_ρ, 2.5)
        set!(Ωfield_0, 2.5)

        aux = biogeochemical_auxiliary_fields(model_ρ.biogeochemistry)

        μP = CUDA.@allowscalar NPDM.PlanktonModels.phytoplankton_growth(1, 1, 1, grid, pl, bgc, model_ρ.tracers, aux)
        Gp = CUDA.@allowscalar NPDM.DetritusModels.grazing(1, 1, 1, grid, Val(:P), pl, bgc, model_ρ.tracers, aux)
        P  = CUDA.@allowscalar model_ρ.tracers.P[1, 1, 1]

        R = pl.carbon_ratio
        γ = pl.phytoplankton_exudation_fraction
        m = pl.phytoplankton_mortality_rate
        η = pl.zooplankton_calcium_carbonate_dissolution

        # ∂CaCO₃ = ρR[(1−η)Gₚ + m P²]   
        expected_CaCO₃ = ρ * R * ((1 - η) * Gp + m * P^2)
        tCaCO₃ = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃), model_ρ.clock, model_ρ.tracers, aux)
        @test tCaCO₃ ≈ expected_CaCO₃

        # DIC calcium carbonate term (isolated against ρ=0): −ρR[(1−γ)μP − η Gₚ]
        tDIC_calcium_carbonate = CUDA.@allowscalar (bgc(1, 1, 1, grid, Val(:DIC), model_ρ.clock, model_ρ.tracers, aux)
                                        - bgc0(1, 1, 1, grid, Val(:DIC), model_0.clock, model_0.tracers, aux))
        tAlk_calcium_carbonate = CUDA.@allowscalar (bgc(1, 1, 1, grid, Val(:Alk), model_ρ.clock, model_ρ.tracers, aux)
                                        - bgc0(1, 1, 1, grid, Val(:Alk), model_0.clock, model_0.tracers, aux))

        expected_DIC_calcium_carbonate = -ρ * R * ((1 - γ) * μP - η * Gp)
        @test tDIC_calcium_carbonate ≈ expected_DIC_calcium_carbonate
        @test tAlk_calcium_carbonate ≈ 2 * expected_DIC_calcium_carbonate

        # tendency-level carbon closure: ρR·∂P + ∂CaCO₃ + (DIC calcium carbonate term) = 0
        ∂P = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:P), model_ρ.clock, model_ρ.tracers, aux)
        @test ρ * R * ∂P + tCaCO₃ + tDIC_calcium_carbonate ≈ 0 atol = 1e-16
    end

    @testset "carbon conservation (explicit CaCO₃ pool)" begin
        ic = ExplicitCalciumCarbonate(grid)
        biogeochemistry = NutrientsPlanktonDetritus(grid;
            plankton = PhytoZoo(grid),
            nutrients = Nutrients(; nitrogen = OceanBioME.N),
            detritus = DissolvedParticulate(grid),
            inorganic_carbon = ic,
            light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100)))

        model = NonhydrostaticModel(grid;
            biogeochemistry, advection = nothing,
            auxiliary_fields = (T = ConstantField(15.0), S = ConstantField(35.0)))

        set!(model, N = 5, P = 1, Z = 0.5, DOM = 1, sPOM = 1, bPOM = 1, DIC = 2000, Alk = 2300, CaCO₃ = 3)

        # explicit element sum: total carbon = DIC + CaCO₃ + organic-pool carbon + living calcium carbonate (ρRP,
        # the calcium carbonate carried by phytoplankton under formation-at-production, weight R(1+ρ) on P).
        bgcu = model.biogeochemistry.underlying_biogeochemistry
        R = OceanBioME.Models.NutrientsPlanktonDetritusModels.carbon_ratio(bgcu.plankton, bgcu)
        ρ = bgcu.plankton.rain_ratio
        total_carbon(m) = CUDA.@allowscalar (m.tracers.DIC[1,1,1] + m.tracers.CaCO₃[1,1,1]
            + R * (1 + ρ) * m.tracers.P[1,1,1]
            + R * (m.tracers.Z[1,1,1] + m.tracers.DOM[1,1,1] + m.tracers.sPOM[1,1,1] + m.tracers.bPOM[1,1,1]))

        C₀ = total_carbon(model)
        initial = get_conservation_values(model.biogeochemistry, model.tracers)

        for _ in 1:100
            time_step!(model, 100)
        end

        C₁ = total_carbon(model)
        final = get_conservation_values(model.biogeochemistry, model.tracers)

        # total carbon is conserved, both by an explicit element sum and via conserved_tracers
        @test isapprox(C₁, C₀; rtol = 1e-8)
        @test isapprox(final.carbon, initial.carbon; rtol = 1e-8)
    end

    @testset "replicates are independent" begin
        model = explicit_calcium_carbonate_model(grid; replicates = 2,
            plankton = PhytoZoo(grid),
            nutrients = Nutrients(; nitrogen = OceanBioME.N),
            detritus = DissolvedParticulate(grid))

        bgc = model.biogeochemistry.underlying_biogeochemistry
        Ω1, Ω2 = bgc.inorganic_carbon.calcium_carbonate_saturation
        aux = biogeochemical_auxiliary_fields(model.biogeochemistry)

        # identical state in both realisations ⇒ identical CaCO₃ tendency. Ω is set by hand
        # (undersaturated ⇒ dissolution active) *after* set!(model, …) since that recomputes it.
        set!(model, N = 5, P = 1, Z = 0.5, DOM = 1, sPOM = 1, bPOM = 1,
             DIC1 = 2000, Alk1 = 2300, CaCO₃1 = 3, DIC2 = 2000, Alk2 = 2300, CaCO₃2 = 3)
        set!(Ω1, 0.5)
        set!(Ω2, 0.5)

        t1 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃1), model.clock, model.tracers, aux)
        t2 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃2), model.clock, model.tracers, aux)
        @test t1 == t2

        # differing CaCO₃ ⇒ differing tendency (dissolution scales with the realisation's own CaCO₃)
        set!(model, CaCO₃2 = 30)
        set!(Ω1, 0.5)
        set!(Ω2, 0.5)
        t1b = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃1), model.clock, model.tracers, aux)
        t2b = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃2), model.clock, model.tracers, aux)
        @test t1b != t2b
    end

    @testset "per-replicate rate parameters (machine precision)" begin
        k1, k2 = 0.2 / day, 0.05 / day    # dissolution rates
        p1, p2 = 0.0, 0.1 / day           # precipitation rates
        m1, m2 = 1.0, 2.0                 # dissolution exponents

        model = explicit_calcium_carbonate_model(grid; replicates = 2,
            calcium_carbonate_dissolution_rate = (k1, k2),
            calcium_carbonate_dissolution_exponent = (m1, m2),
            calcium_carbonate_precipitation_rate = (p1, p2),
            calcium_carbonate_precipitation_exponent = 2)

        bgc = model.biogeochemistry.underlying_biogeochemistry
        ic  = bgc.inorganic_carbon

        # scalars broadcast, tuples are kept per replicate and in order
        @test ic.calcium_carbonate_dissolution_rate       == (k1, k2)
        @test ic.calcium_carbonate_dissolution_exponent   == (m1, m2)
        @test ic.calcium_carbonate_precipitation_rate     == (p1, p2)
        @test ic.calcium_carbonate_precipitation_exponent == (2.0, 2.0)

        Ω1, Ω2 = ic.calcium_carbonate_saturation
        aux = biogeochemical_auxiliary_fields(model.biogeochemistry)

        # identical state, undersaturated ⇒ only the dissolution parameters differ
        CaCO₃, Ωval = 3.0, 0.4
        set!(model, DIC1 = 2000, Alk1 = 2300, CaCO₃1 = CaCO₃,
                    DIC2 = 2000, Alk2 = 2300, CaCO₃2 = CaCO₃)
        set!(Ω1, Ωval)
        set!(Ω2, Ωval)

        t1 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃1), model.clock, model.tracers, aux)
        t2 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃2), model.clock, model.tracers, aux)

        # each replicate must use *its own* k and m (Abiotic plankton ⇒ no biological production)
        @test t1 ≈ -k1 * (1 - Ωval)^m1 * CaCO₃
        @test t2 ≈ -k2 * (1 - Ωval)^m2 * CaCO₃
        @test t1 != t2

        # DIC/Alk pick up the same per-replicate rates, 1x and 2x
        d1 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:DIC1), model.clock, model.tracers, aux)
        d2 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:DIC2), model.clock, model.tracers, aux)
        a2 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:Alk2), model.clock, model.tracers, aux)
        @test d1 ≈  k1 * (1 - Ωval)^m1 * CaCO₃
        @test d2 ≈  k2 * (1 - Ωval)^m2 * CaCO₃
        @test a2 ≈ 2k2 * (1 - Ωval)^m2 * CaCO₃

        # supersaturated ⇒ only the precipitation parameters differ, and p1 = 0 disables it
        set!(Ω1, 2.5)
        set!(Ω2, 2.5)
        s1 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃1), model.clock, model.tracers, aux)
        s2 = CUDA.@allowscalar bgc(1, 1, 1, grid, Val(:CaCO₃2), model.clock, model.tracers, aux)
        @test s1 == 0
        @test s2 ≈ p2 * (2.5 - 1)^2 * CaCO₃

        # a wrong-length tuple is rejected up front rather than in a kernel
        @test_throws ArgumentError ExplicitCalciumCarbonate(grid; replicates = 2,
            calcium_carbonate_dissolution_rate = (k1, k2, k1))
    end
end
