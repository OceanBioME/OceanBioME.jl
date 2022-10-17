    using Test
using OceanBioME: PISCES
using Oceananigans
using Oceananigans.Fields: TracerFields, fill_halo_regions!

#initial values
P = 4e-8; Z = 4e-8; D = 4e-8; M = 4e-8; Chlᴾ = 6e-10; Chlᴰ = 6e-10; Feᴾ = 2.5e-7; Feᴰ = 2.5e-7; Siᴰ = 5.5e-9; NO₃ = 11.4; NH₄ = 0.05; PO₄ = 11.4; Fe = 11.4; Si =11.4; T=15; DOC = 0.0; POC = 0.0; GOC = 0.0; Feᴾᴼ = 0.0; Feᴳᴼ = 0.0; Siᴾ = 0.0; CaCO₃ = 0.0; DIC = 2200; Alk = 2400; O₂ = 240; Si̅ = 12
PAR¹ = 10; PAR² = 10; PAR³ = 10; zₑᵤ = 5; zₘₓₗ = 5; D_dust = 1; S=35

L_day(t) = 0.5
ϕ = 10
params = merge(PISCES.defaults, (; ϕ, L_day)) 

#need this for discrete forcing, should return same position etc. as inputted in continuos functions
clock = Clock(time=0)
grid = RectilinearGrid(size=(1, 1, 1), extent=(2, 2, 20), topology=(Periodic, Periodic, Periodic))
fields = (;P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₑᵤ, zₘₓₗ, Si̅, D_dust)
model_fields = TracerFields(keys(fields), grid)
for field in keys(fields)
    model_fields[field] .= @eval $field
end
fill_halo_regions!(model_fields)
@testset "Realistic Values" begin
    for tracer in PISCES.tracers
        @info "Testing $tracer"
        fname = Symbol(tracer, "_forcing")
        if PISCES.discrete_forcing[tracer]
            d = @eval PISCES.$fname(1, 1, 1, grid, clock, model_fields, params)
        else
            d = @eval PISCES.$fname(1, 1, -10, 0, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₑᵤ, zₘₓₗ, Si̅, D_dust, params)
        end
        @test !isnan(d)
    end
end

#initial values
P = 0.0; Z = 0.0; D = 0.0; M = 0.0; Chlᴾ = 0.0; Chlᴰ = 0.0; Feᴾ = 0.0; Feᴰ = 0.0; Siᴰ = 0.0; NO₃ = 0.0; NH₄ = 0.0; PO₄ = 0.0; Fe = 0.0; Si =0.0; T=15; DOC = 0.0; POC = 0.0; GOC = 0.0; Feᴾᴼ = 0.0; Feᴳᴼ = 0.0; Siᴾ = 0.0; CaCO₃ = 0.0; DIC = 0.0; Alk = 0.0; O₂ = 0.0; Si̅ = 0.0
PAR¹ = 0.0; PAR² = 0.0; PAR³ = 0.0; zₑᵤ = 5; zₘₓₗ = 5; D_dust = 0.0; S=35

fields = (;P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₑᵤ, zₘₓₗ, Si̅, D_dust)
model_fields = TracerFields(keys(fields), grid)
for field in keys(fields)
    model_fields[field] .= @eval $field
end
fill_halo_regions!(model_fields)

@testset "Zeros" begin
    for tracer in PISCES.tracers
        @info "Testing $tracer"
        fname = Symbol(tracer, "_forcing")
        if PISCES.discrete_forcing[tracer]
            d = @eval PISCES.$fname(1, 1, 1, grid, clock, model_fields, params)
        else
            d = @eval PISCES.$fname(1, 1, -10, 0, P, D, Chlᴾ, Chlᴰ, Feᴾ, Feᴰ, Siᴰ, Z, M, DOC, POC, GOC, Feᴾᴼ, Feᴳᴼ, Siᴾ, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, PAR¹, PAR², PAR³, T, S, zₑᵤ, zₘₓₗ, Si̅, D_dust, params)
        end
        @test !isnan(d)
    end
end