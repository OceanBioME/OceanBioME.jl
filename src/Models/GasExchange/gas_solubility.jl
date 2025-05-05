"""
    PartiallySolubleGas(; air_concentration, solubility)

Parameterises the available concentration of a gas dissolving in water in the form ``\\alpha C_a``
where ``\alpha`` is the Ostwald solubility coeffieient and ``C_a`` is the concentration in the air.
"""
struct PartiallySolubleGas{AC, S}
    air_concentration :: AC
           solubility :: S

    function PartiallySolubleGas(FT = Float64;
                                 air_concentration,
                                 solubility :: S) where S
        
        air_concentration = normalise_surface_function(air_concentration; FT)
        
        AC = typeof(air_concentration)

        return new{AC, S}(air_concentration, solubility)
    end
end

@inline surface_value(gs::PartiallySolubleGas, i, j, grid, clock, model_fields) = 
    surface_value(gs.air_concentration, i, j, grid, clock) * surface_value(gs.solubility, i, j, grid, clock, model_fields)

"""
    Wanninkhof92Solubility

Parameterises the Ostwald solubility coefficient as given in Wanninkhof, 1992.
"""
struct Wanninkhof92Solubility{FT}
    A1 :: FT
    A2 :: FT
    A3 :: FT
    B1 :: FT
    B2 :: FT
    B3 :: FT
end

function surface_value(sol::Wanninkhof92Solubility, i, j, grid, clock, model_fields)
    FT = eltype(grid)

    Tk = @inbounds model_fields.T[i, j, grid.Nz] + convert(FT, 273.15)
    S = @inbounds model_fields.S[i, j, grid.Nz]

    Tk_100 = Tk / convert(FT, 100)

    β = exp(sol.A1 + sol.A2 / Tk_100 + sol.A3 * log(Tk_100) + S * (sol.B1 + sol.B2 * Tk_100 + sol.B2 * Tk_100^convert(FT, 2)))

    return β / Tk
end

OxygenSolubility(FT = Float64; A1 = -58.3877, A2 = 85.8079, A3 = 23.8439, B1 = -0.034892, B2 = 0.015578, B3 = -0.0019387) =
    Wanninkhof92Solubility{FT}(A1, A2, A3, B1, B2, B3)