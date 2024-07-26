using OceanBioME, Test, CUDA, Oceananigans, JLD2, Documenter

using OceanBioME.Sediments: SimpleMultiG, InstantRemineralisation
using Oceananigans.Units

using OceanBioME.Sediments: sediment_tracers, sediment_fields
using Oceananigans: Field
using Oceananigans.Fields: TracerFields

using Oceananigans.Operators: volume, Azᶠᶜᶜ

using OceanBioME.LOBSTERModel: VariableRedfieldLobster

using CairoMakie
using CairoMakie: record

Z = FieldTimeSeries("LDNtesting.jld2", "b")
N = FieldTimeSeries("LDNtesting.jld2", "Z")
P = FieldTimeSeries("LDNtesting.jld2", "P")

xc, yc, zc = nodes(Z)

times = Z.times

print(nodes(Z))
9
Z_lims = (findmin(Z.data)[1], findmax(Z.data)[1])
N_lims = (findmin(N.data)[1], findmax(N.data)[1])
P_lims = (findmin(P.data)[1], findmax(P.data)[1])

n = Observable(1)

Zₙ = @lift interior(Z[$n], :, 2, :)
Nₙ = @lift interior(N[$n], :, 2, :)
Pₙ = @lift interior(P[$n], :, 2, :)

println(Pₙ)

fig = Figure(size = (1000, 520), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[2, 3], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 150) #, yticks = [-400, -200, 0])
ax1 = Axis(fig[1, 1]; title = "Zooplankton concentration\n (mmol N / m³)", axis_kwargs...)
ax2 = Axis(fig[1, 3]; title = "NO₃ concentration\n (mmol N / m³)",axis_kwargs...)
ax3 = Axis(fig[1, 5]; title = "Phytoplankton concentration\n (mmol N / m³)", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, Zₙ, colorrange = Z_lims, colormap = Reverse(:bamako), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, Nₙ, colorrange = N_lims, colormap = Reverse(:lajolla), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1)
Colorbar(fig[1, 4], hm2)
Colorbar(fig[1, 6], hm3)

rowgap!(fig.layout, 0)

record(fig, "LDNtesting.gif", 1:length(times)) do i
    n[] = i
end 