using Test
using OceanBioME: Boundaries

@testset "Gas exchange" begin
    # approximatly correct sized pCO₂ from DIC and ALK
    pCO₂ = Boundaries.pCO₂(2220, 2500, 15+273.15, 35, Boundaries.defaults.airseaflux)
    @test pCO₂ ≈ 400 atol = 100

    # correct sized flux based on pCO₂ = 355 atm and NPP of ∼3mg C/m³/day in top 100m (from Copernicus data) -> ∼0.0003 mmol C/m²/s
    flux = Boundaries.airseaflux(0, 0, 0, 15.0, 35.0, 350.0, merge(Boundaries.defaults.airseaflux, (gas=:CO₂, )))
    @test flux ≈ -0.0003 atol = 0.0001

    # Oxygen exchange (since most gases use a slightly different formulation), must also be same order as NPP*Rd_resp ≈ NPP/0.8 -> ∼0.0004 mmol O₂/m²/s
    flux = Boundaries.airseaflux(0, 0, 0, 15.0, 35.0, 270.0, merge(Boundaries.defaults.airseaflux, (gas=:O₂, )))
    @test flux ≈ -0.0004 atol = 0.0001
end