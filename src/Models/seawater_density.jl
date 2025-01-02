using GibbsSeaWater: gsw_rho, gsw_sa_from_sp
using SeawaterPolynomials.TEOS10: _ρ, τ, s, ζ

"""
    teos10_density(T, Sp, Pbar = 0, lon = 0, lat = 0)

Returns sea water density computed by `GibbsSeaWater.jl` from
`T` in degrees centigrade, `Sp` in PSU,  `P`ressure in bar. Location
`lat`itude and `lon`gitude are required for accurate conversion from
practical to absolute salinity.

This function *will not* run on GPU so can not be used as the default.
"""
@inline function teos10_density(T, Sp, Pbar = 0, lon = 0, lat = 0)
    Pdbar = 10 * Pbar

    Sa = gsw_sa_from_sp(Sp, Pdbar, lon, lat)

    return gsw_rho(Sa, T, Pdbar)
end

"""
    teos10_polynomial_approximation(T, Sp, Pbar = 0, lon = 0, lat = 0)

Returns sea water density computed by `SeawaterPolynomials.jl` TEOS10
polynomial approximation from `T` in degrees centigrade, `S` in PSU,
`P`ressure in bar. This approximation also requires us to approximate
geopotential depth from pressure at 1 m / dbar, and absolute salinity
as practial salinity.

This method is less accurate than `teos10_density` but will run on GPU.
"""
@inline function teos10_polynomial_approximation(T, Sp, Pbar = 0, lon = 0, lat = 0)
    Z = 10 * Pbar # 10 m / bar

    Sa = Sp # error is *very* small in the pH and pCO₂ from this approximation

    return _ρ(τ(T), s(Sa), ζ(Z))
end