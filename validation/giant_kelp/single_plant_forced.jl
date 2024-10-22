using CairoMakie, Oceananigans.Units, OceanBioME, NetCDF, JLD2, Statistics, Interpolations

no3 = mean(ncread("../GiantKelpGrowth/cmems_nitrate.nc", "no3")[1, 1, :, :], dims=1)[1, :]
no3_times = ncread("../GiantKelpGrowth/cmems_nitrate.nc", "time")

no3_times .-= no3_times[1] # now seconds since 01/01/1996

temp = mean(ncread("../GiantKelpGrowth/cmems_phys.nc", "thetao")[1, 1, :, :], dims=1)[1, :]
#sal = ncread("cmems_phys.nc", "so")
phy_times = ncread("../GiantKelpGrowth/cmems_phys.nc", "time") # times go 01/01/1996 to 31/12/2018

phy_times .-= phy_times[1] # now seconds since 01/01/1996

@load "../GiantKelpGrowth/par_ts.jld2" PAR_itp # 1998/01/01 to 2007/12/31

# cut no3 and temp to 1998/01/01 to 2007/12/31

no3 = no3[732:4382]
no3_times = no3_times[732:4382]

no3_times .-= no3_times[1]

temp = temp[732:4382]
phy_times = phy_times[732:4382]

phy_times .-= phy_times[1]

no3_itp = LinearInterpolation(no3_times, no3)
temp_itp = LinearInterpolation(phy_times, temp)

# make the model
model = GiantKelp()

A, N, C = [1.0], [0.02], [0.01]

t = [366Units.days] .+ 240Units.days .+ 366Units.days

te = t[]+60Units.days#3650Units.days
Δt = 10Units.minutes

nt = floor(Int, (te-t[]) / Δt)

A_ts, N_ts, C_ts = zeros(nt), zeros(nt), zeros(nt)

times = [t[]+Δt:Δt:te;] .- t[] .- Δt

for n in 1:nt
    A_ts[n] = A[]
    N_ts[n] = N[]
    C_ts[n] = C[]

    t[] += Δt

    T = temp_itp(t[])
    NO₃ = no3_itp(t[])
    NH₄ = 0.29 * NO₃
    PAR = PAR_itp(t[])# * 10^6 / day# mol / m^2 / day to μ mol / m^2 / s #(3.99e-10 * 545e12) / day * exp(-2*.2) # mol / m^2 / day -> w / m² / s
    
    A[] += model(Val(:A), t[], A[], N[], C[], 1, 0, 0, T, NO₃, NH₄, PAR) * Δt
    N[] += model(Val(:N), t[], A[], N[], C[], 1, 0, 0, T, NO₃, NH₄, PAR) * Δt
    C[] += model(Val(:C), t[], A[], N[], C[], 1, 0, 0, T, NO₃, NH₄, PAR) * Δt

    C[] < 0.01 && @error "C went below min"
end

