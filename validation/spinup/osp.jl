using OceanBioME, Oceananigans, Statistics
using NetCDF, JLD2, Interpolations, Oceananigans.Units
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Simulations: UpdateStateCallsite

# settings
path = "validation/spinup/"

# load data
PAR = ncread(path*"PAR2.nc", "PAR")[1:30] ./ day ./ 4.56 .* 1e6
PAR_times = ncread(path*"PAR2.nc", "times")[1:30]

PAR_itp = linear_interpolation(PAR_times, PAR; extrapolation_bc = Interpolations.Periodic())
const mean_PAR = mean(PAR)
@inline PAR_static(t) = mean_PAR

@kwdef struct SolarShortwave{MI, FT, IT} <: Function
     mean_intensity :: MI = PAR_static#PAR_itp
     local_meridiem :: FT = 0hour
           latitude :: FT = 50.1
 initial_day_number :: IT = 188
end

adapt_structure(to, solar::SolarShortwave) = solar

@inline (ss::SolarShortwave)(args...) = @inbounds ss(args[end])

@inline function (solar_shortwave::SolarShortwave)(t_)
    tₘ = solar_shortwave.local_meridiem
    L = solar_shortwave.latitude
    N₀ = solar_shortwave.initial_day_number
    Ī = solar_shortwave.mean_intensity(t_)

    t = t_ - tₘ - 12hours
    N = N₀ + floor(Int, t / 24hours)

    δ = 23.45 * sind(360 * (284 + N) / 365)

    # tells us the day length
    a = sind(L) * sind(δ)

    # normalise the mean intensity to get the peak intensity given the day length
    I₀ = Ī * (1 + a) / (sin(acos(-a)))
    
    # reflectivity
    # from https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.4319/lo.1990.35.8.1657
    # maybe change to wind speed dependant formula
    h = t * 15 / hours
    α = 90 - asind(sind(L) * sind(δ) + cosd(L) * cosd(δ) * cosd(h))
    αᵣ = asind(sind(α) / 1.431)
    ρ = 1/2 * ((sind(α-αᵣ)/sind(α+αᵣ))^2 / 2 + (tand(α-αᵣ)/tand(α+αᵣ))^2)

    return max(0, (1 - ρ)) * I₀ * max(0, (a + cos(2π * t / 24hours)) / (1 + a))
end

bgc_z = reverse(-ncread(path*"P_0707_2709.nc", "depth"))
P = linear_interpolation(bgc_z, reverse(ncread(path*"P_0707_2709.nc", "phyc")[1, 1, :, 1] ./ 6.56))
Z = linear_interpolation(bgc_z, reverse(ncread(path*"Z_0707_2709.nc", "zooc")[1, 1, :, 1] ./ 6.56))
NO₃ = linear_interpolation(bgc_z, reverse(ncread(path*"NO3_0707_2709.nc", "no3")[1, 1, :, 1]))
Fe = linear_interpolation(bgc_z, reverse(ncread(path*"Fe_0707_2309.nc", "fe")[1, 1, :, 1]))
DIC = linear_interpolation(bgc_z, reverse(ncread(path*"chem_0707_2709.nc", "dissic")[1, 1, :, 1] .* 1000))
Alk = linear_interpolation(bgc_z, reverse(ncread(path*"chem_0707_2709.nc", "talk")[1, 1, :, 1] .* 1000))

# these are from 2024
T = linear_interpolation(bgc_z, reverse(ncread(path*"temp_0707_2709.nc", "thetao")[1, 1, :, 1]))
S = linear_interpolation(bgc_z, reverse(ncread(path*"salinity_0707_2709.nc", "so")[1, 1, :, 1]))

# model setup
grid = RectilinearGrid(topology = (Oceananigans.Flat, Oceananigans.Flat, Oceananigans.Bounded), size = (128, ), extent = (500, ))

biogeochemistry = LOBSTER(; grid,
                            surface_photosynthetically_active_radiation = PAR_itp,
                            carbonates = true,
                            iron_half_saturation = 2e-4#2.1e-4#2.5e-4,#5e-4,#1e-5, # MARBL #3e-6, # PISCES code namelists #0.001, # PISCES
                            )

nudging_timescale = 15days

NO₃_restoring = Forcing((z, t, NO₃, NO₃_initial) -> (NO₃_initial(z) - NO₃) / nudging_timescale,#min(2days, ifelse(t>20days, (t-20days)/80days * .9days + .1days, .1days)), 
                        parameters = NO₃, field_dependencies = (:NO₃, ))

Fe_restoring = Forcing((z, t, Fe, Fe_initial) -> (Fe_initial(z) - Fe) / nudging_timescale,#/ min(2days, ifelse(t>20days, (t-20days)/80days * .9days + .1days, .1days)), 
                        parameters = Fe, field_dependencies = (:Fe, ))

CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition(;air_concentration = 414.0)

DIC_bcs = FieldBoundaryConditions(top = CO₂_flux)

Fe_deposition = FluxBoundaryCondition(6e-12 / 55.845 * 1000)#6pg/m^2/s to mmol Fe / m^2 / s (= mmol Fe / m^3 * m/s)

Fe_bcs = FieldBoundaryConditions(top = Fe_deposition)

N_flux_out = Ref(0.0)
N_flux_in = Ref(0.0)

Fe_flux_in = Ref(0.0)

model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              tracers = (:T, :S),
                              forcing = (; NO₃ = NO₃_restoring, Fe = Fe_restoring),
                              boundary_conditions = (; DIC = DIC_bcs, Fe_bcs),
                              closure = ScalarDiffusivity(κ = (z, t) -> 1e-5 + 1e-2 * (1 - min(1, max(0, -(z+6)/5)))))

# simulation
set!(model, P = (z)->P(z), Z = (z)->Z(z), NO₃ = (z)->NO₃(z), DIC = (z)->DIC(z), Alk = (z)->Alk(z), T= (z)->T(z) , S = (z)->S(z), Fe = (z)->Fe(z))

N₀ = sum(model.tracers.P) + sum(model.tracers.Z) + sum(model.tracers.NO₃)

simulation = Simulation(model, Δt = 3minutes, stop_time = 500days)

prog(sim) = @info prettytime(sim)*" in "*prettytime(sim.run_wall_time)

add_callback!(simulation, prog, IterationInterval(1000))

simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers; schedule = TimeInterval(1day), filename = "spinup_osp", overwrite_existing = true)

N_out_record = Float64[]
N_in_record = Float64[]
Fe_in_record = Float64[]

function update_flux(sim)
    model = sim.model
    # I think this is maybe not technically correct or maybe just not the most correct way to compute this
    @inbounds begin
        sPOM_loss = model.tracers.sPOM[1, 1, 1] * biogeochemical_drift_velocity(model.biogeochemistry, Val(:sPOM)).w[1, 1, 1]
        bPOM_loss = model.tracers.bPOM[1, 1, 1] * biogeochemical_drift_velocity(model.biogeochemistry, Val(:bPOM)).w[1, 1, 1]
# mmol N / m^2 / s
        N_flux_out[] = -(sPOM_loss + bPOM_loss)

        N_flux_in[] = sum([500 / 128 * (NO₃(z) - model.tracers.NO₃[1, 1, k]) / nudging_timescale for (k, z) in enumerate(znodes(model.tracers.NO₃))])
        Fe_flux_in[] = sum([500 / 128 * (NO₃(z) - model.tracers.Fe[1, 1, k]) / nudging_timescale for (k, z) in enumerate(znodes(model.tracers.Fe))])

        push!(N_out_record, N_flux_out[])
        push!(N_in_record, N_flux_in[])
        push!(Fe_in_record, Fe_flux_in[])
    end

    return nothing
end

add_callback!(simulation, update_flux, IterationInterval(10))


run!(simulation)

using CairoMakie

fds = FieldDataset("spinup_osp.jld2")

fig = Figure();
ax = Axis(fig[1, 1])
n = Observable(1)

fds[:NO₃] ./= 10
fds[:Fe] .*= 1000

for tracer in (:P, :Z, :DOM, :sPOM, :bPOM, :NO₃, :NH₄, :Fe)
    lines!(ax, (@lift fds[tracer][$n]), label = "$tracer"*ifelse(tracer==:NO₃, " (÷10)", ifelse(tracer==:Fe, " (x1000)", "")))
end

xlims!(ax, 0, 3)

axislegend(ax, position = :rb)

record(fig, "spinup_static_even_high_kFe.mp4", 1:length(fds[:P].times); framerate = 10) do i; 
    n[] = i
end

fig2 = Figure()

ax2 = Axis(fig2[1, 1])

lines!(ax2, fds[:P][1], color = Makie.wong_colors()[1], linestyle = :dash, label = "P (mmol N /m³)")
lines!(ax2, fds[:P][end], color = Makie.wong_colors()[1])

lines!(ax2, fds[:NO₃][1], color = Makie.wong_colors()[2], linestyle = :dash, label = "NO₃ (÷10 mmol N/m³)")
lines!(ax2, fds[:NO₃][end], color = Makie.wong_colors()[2])

lines!(ax2, fds[:Fe][1], color = Makie.wong_colors()[3], linestyle = :dash, label = "Fe (μmol Fe/m³)")
lines!(ax2, fds[:Fe][end], color = Makie.wong_colors()[3])

# TODO: compare to chlorophyll obs