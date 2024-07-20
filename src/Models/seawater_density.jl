using GibbsSeaWater: gsw_rho, gsw_sa_from_sp

"""
    seawater_density(T, S, P=0)

Returns sea water density computed by `GibbsSeaWater.jl` from
`T` in degrees centigrade, `S` in PSU, and pressure in decibar.
"""
@inline function seawater_density(T, Sp, P = 0, lon = 0, lat = 0)
    Sa = gsw_sa_from_sp(Sp, P, lon, lat)
    
    return gsw_rho(Sa, T, P)
end