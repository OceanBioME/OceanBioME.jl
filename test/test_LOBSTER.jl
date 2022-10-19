using Test
using OceanBioME: LOBSTER
using Oceananigans
using Oceananigans.Fields: TracerFields, fill_halo_regions!

#initial values
P = 0.004*rand(); Z = 0.004*rand(); D = 0.002*rand(); DD = 0.002*rand(); NO₃ = 11.4*rand(); NH₄ = 0.05*rand(); DOM = 0.0*rand(); DIC = 2200*rand(); ALK = 2400*rand(); OXY = 240*rand(); 
PAR = 10*rand(); 

function test_forcing(val)
    @test !isnan(val)
    return val
end

clock = Clock(time=0)
grid = RectilinearGrid(size=(1, 1, 1), extent=(2, 2, 2), topology=(Periodic, Periodic, Periodic))
fields = (;P, Z, D, DD, NO₃, NH₄, DOM, PAR)
model_fields = TracerFields(keys(fields), grid)
for field in keys(fields)
    model_fields[field] .= @eval $field
end
fill_halo_regions!(model_fields)

@testset "Realistic Values" begin
    dP = test_forcing(LOBSTER.P_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(1, 1, 1, grid, clock, model_fields, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(1, 1, 1, grid, clock, model_fields, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR, LOBSTER.defaults))

    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0 atol = max([eps(x) for x in [P, Z, D, DD, DOM, NO₃, NH₄]]...)
end

P = 0.0; Z = 0.0; D = 0.0; DD = 0.0; NO₃ = 0.0; NH₄ = 0.0; DOM = 0.0; DIC = 0.0; ALK = 0.0; OXY = 0.0; 
PAR = 0.0; 

fields = (;P, Z, D, DD, NO₃, NH₄, DOM, PAR)
model_fields = TracerFields(keys(fields), grid)
for field in keys(fields)
    model_fields[field] .= @eval $field
end
fill_halo_regions!(model_fields)
@testset "Zeros" begin
    dP = test_forcing(LOBSTER.P_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(1, 1, 1, grid, clock, model_fields, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(1, 1, 1, grid, clock, model_fields, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, LOBSTER.defaults))

    println("$dP + $dZ + $dD + $dDD + $dDOM + $dNO₃ + $dNH₄")
    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0 atol = max([eps(x) for x in [P, Z, D, DD, DOM, NO₃, NH₄]]...)
end
