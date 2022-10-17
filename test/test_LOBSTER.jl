using Test
using OceanBioME: LOBSTER
using Oceananigans
using Oceananigans.Fields: TracerFields, fill_halo_regions!

#initial values
P = 0.004; Z = 0.004; D = 0.002; DD = 0.002; NO₃ = 11.4; NH₄ = 0.05; DOM = 0.0; DIC = 2200; ALK = 2400; OXY = 240; 
PAR = 10; 

function test_forcing(val)
    @test !isnan(val)
    return val
end

@testset "Realistic Values" begin
    dP = test_forcing(LOBSTER.P_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, PAR, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR, LOBSTER.defaults))

    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0 atol = max([eps(x) for x in [P, Z, D, DD, DOM, NO₃, NH₄]]...)
end

@testset "Zeros" begin
    dP = test_forcing(LOBSTER.P_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, LOBSTER.defaults))

    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0
end
