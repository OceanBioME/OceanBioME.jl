using Test
using OceanBioME: LOBSTER
using Oceananigans
using Oceananigans.Fields: TracerFields, fill_halo_regions!

#initial values
P = 0.004*rand(); Z = 0.004*rand(); D = 0.002*rand(); DD = 0.002*rand(); NO₃ = 11.4*rand(); NH₄ = 0.05*rand(); DOM = 0.0*rand(); DIC = 2200*rand(); ALK = 2400*rand(); OXY = 240*rand(); 
Dᶜ = D*LOBSTER.defaults.Rd_phy; DDᶜ = DD*LOBSTER.defaults.Rd_phy;
PAR = 10*rand(); 

function test_forcing(val)
    @test !isnan(val)
    return val
end

@testset "Realistic Values" begin
    dP = test_forcing(LOBSTER.P_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dDᶜ = test_forcing(LOBSTER.Dᶜ_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dDDᶜ = test_forcing(LOBSTER.DDᶜ_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0, 0, -1, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0, 0, 0, 0, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR, LOBSTER.defaults))

    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0 atol = max([eps(x) for x in [P, Z, D, DD, DOM, NO₃, NH₄]]...)
    @test dDᶜ ≈ dD*LOBSTER.defaults.Rd_phy
    @test dDDᶜ ≈ dDD*LOBSTER.defaults.Rd_phy
end

@testset "Zeros" begin
    dP = test_forcing(LOBSTER.P_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dZ = test_forcing(LOBSTER.Z_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dD = test_forcing(LOBSTER.D_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dDD = test_forcing(LOBSTER.DD_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))    
    dDᶜ = test_forcing(LOBSTER.D_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dDDᶜ = test_forcing(LOBSTER.DD_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dDOM = test_forcing(LOBSTER.DOM_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dNO₃ = test_forcing(LOBSTER.NO₃_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dNH₄ = test_forcing(LOBSTER.NH₄_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))

    dDIC = test_forcing(LOBSTER.DIC_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    dALK = test_forcing(LOBSTER.ALK_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))

    dOXY = test_forcing(LOBSTER.OXY_forcing(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, LOBSTER.defaults))
    
    @test (dP + dZ + dD + dDD + dDOM + dNO₃ + dNH₄) ≈ 0 atol = max([eps(x) for x in [P, Z, D, DD, DOM, NO₃, NH₄]]...)
    @test dDᶜ ≈ dD*LOBSTER.defaults.Rd_phy
    @test dDDᶜ ≈ dDD*LOBSTER.defaults.Rd_phy
end
