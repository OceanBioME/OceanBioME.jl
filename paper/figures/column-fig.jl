using JLD2, CairoMakie, Statistics
using Oceananigans.Units

@inline t_function(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

# could change this to FieldTimeSeries but would be harder to stick the arrays together (probably could also write to the same file but too late)
P, NO₃, DIC, Alk, sPOC, bPOC = ntuple(n -> zeros(50, 4 * 365), 6)

times = ones(4 * 365) .* NaN

file = jldopen("column.jld2")
iterations = keys(file["timeseries/t"])

for (idx, it) in enumerate(iterations)
    P[:, idx] = file["timeseries/P/$it"][1, 1, 1:50]
    NO₃[:, idx] = file["timeseries/NO₃/$it"][1, 1, 1:50]
    DIC[:, idx] = file["timeseries/DIC/$it"][1, 1, 1:50]
    Alk[:, idx] = file["timeseries/Alk/$it"][1, 1, 1:50]
    sPOC[:, idx] = file["timeseries/sPOC/$it"][1, 1, 1:50]
    bPOC[:, idx] = file["timeseries/bPOC/$it"][1, 1, 1:50]
    times[idx] = file["timeseries/t/$it"]
end

close(file)

initial_offset = length(iterations) - 1

file = jldopen("kelp.jld2")
iterations = keys(file["timeseries/t"])[2:end]

for (idx, it) in enumerate(iterations)
    P[:, idx + initial_offset] = file["timeseries/P/$it"][1, 1, 1:50]
    NO₃[:, idx + initial_offset] = file["timeseries/NO₃/$it"][1, 1, 1:50]
    DIC[:, idx + initial_offset] = file["timeseries/DIC/$it"][1, 1, 1:50]
    Alk[:, idx + initial_offset] = file["timeseries/Alk/$it"][1, 1, 1:50]
    sPOC[:, idx + initial_offset] = file["timeseries/sPOC/$it"][1, 1, 1:50]
    bPOC[:, idx + initial_offset] = file["timeseries/bPOC/$it"][1, 1, 1:50]
    times[idx + initial_offset] = file["timeseries/t/$it"]
end

times .-= times[1]

close(file)

file = jldopen("kelp_particles.jld2")
iterations = keys(file["timeseries/t"])[2:end-183]

A, N, C = ntuple(n -> ones(5, length(iterations)) .* NaN, 3)
kelp_times = zeros(length(iterations))

for (idx, it) in enumerate(iterations)
    particles = file["timeseries/particles/$it"]
    A[:, idx] = particles.A
    N[:, idx] = particles.N
    C[:, idx] = particles.C
    kelp_times[idx] = file["timeseries/t/$it"]
end

kelp_times .-= 3years

close(file)

air_sea_CO₂_flux = similar(times)
carbon_export = similar(times)
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[50, i], Alk[50, i], t_function(1, 1, 0, t), 35)
    carbon_export[i] = (sPOC[end - 20, i] * 3.47e-5 + bPOC[end - 20, i] * 200/day) 
end

zs = [-198:4:-2;]

fig = Figure(resolution = (1200, 1000))

axP = Axis(fig[1:3, 1:2], xlabel = "Time (years)", ylabel = "Depth (m)", title = "(a) Phytoplankton concentration (mmol N / m³)", limits = (0, 3, -160, 0))

hmP = heatmap!(axP, times ./ year, zs[10:50], P[10:50, :]', colormap = Reverse(:batlow))
Colorbar(fig[1:3, 3], hmP)

lines!(axP, [2 - 30/365, 2 - 30/365], [zs[9], 0], color=:black)

axN = Axis(fig[4:6, 1:2], xlabel = "Time (years)", ylabel = "Depth (m)", title = "(b) Nutrient concentration (mmol N / m³)", limits = (0, 3, -160, 0))

hmN = heatmap!(axN, times ./ year, zs[10:50], NO₃[10:50, :]', colormap = Reverse(:batlow))
Colorbar(fig[4:6, 3], hmN)

lines!(axN, [2 - 30/365, 2 - 30/365], [zs[9], 0], color=:black)

axS = Axis(fig[7:8, 1:4], xlabel = "Time (years)", ylabel = "Carbon Flux (kg CO₂ / m² / year)", limits = (0, times[end - 2] /years - 1, 1.1 * (12 + 16 * 2) * year /(1000 * 1000) * min(minimum(air_sea_CO₂_flux[isfinite.(air_sea_CO₂_flux)]), minimum(-carbon_export[isfinite.(carbon_export)])), 1.2 * (12 + 16 * 2) * year /(1000 * 1000) * max(maximum(air_sea_CO₂_flux[isfinite.(air_sea_CO₂_flux)]), maximum(-carbon_export[isfinite.(carbon_export)]))))

lnAS = lines!(axS, times ./ year, air_sea_CO₂_flux .* (12 + 16 * 2) .* year /(1000 * 1000), label = "Air-sea exchange")
lnSI = lines!(axS, times ./ year, - carbon_export .* (12 + 16 * 2) .* year /(1000 * 1000), label = "Sinking export")

Legend(fig[7:8, 1:4], [lnAS, lnSI], ["Air-sea CO₂ exchange", "Sinking export"], halign = :left, valign = :bottom)

lines!(axS, [2 - 30/365, 2 - 30/365], [1, -2], color=:black)

Ā = mean(A, dims=1)[1, :]

axA = Axis(fig[1:2, 4], xlabel = "Time (years)", ylabel = "Frond area (dm² / frond)", title = "(c) Kelp growth", xticks = [2:0.2:3;], limits = (2 - 30/365, 3, minimum(Ā)*0.95, maximum(Ā)*1.05))

lines!(axA, kelp_times ./ year, Ā)

ΣC = sum((C .+ 0.2) .* A .* 0.5 .* 100 * (12 + 16 * 2) / (12 * 1000), dims = 1)[1, :]

axC = Axis(fig[3:4, 4], xlabel = "Time (years)", ylabel = "Carbon stored (kg CO₂ / m²)", xticks = [2:0.2:3;], limits = (2 - 30/365, 3, minimum(ΣC)*0.95, maximum(ΣC)*1.05))

lines!(axC, kelp_times ./ year, ΣC)

ΣN = sum((N .+ 0.0146) .* A .* 0.5 .* 100 / 14, dims = 1)[1, :]

axN = Axis(fig[5:6, 4], xlabel = "Time (years)", ylabel = "Nitrogen stored (mmol N / m²)", xticks = [2:0.2:3;], limits = (2 - 30/365, 3, minimum(ΣN)*0.95, maximum(ΣN)*1.05))

lines!(axN, kelp_times ./ year, ΣN)

save("paper/column_example.png", fig, px_per_unit = 4)