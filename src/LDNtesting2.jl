T = FieldTimeSeries("buoyancy_front.jld2", "T")
N = FieldTimeSeries("buoyancy_front.jld2", "N")
P = FieldTimeSeries("buoyancy_front.jld2", "P")

xc, yc, zc = nodes(T)

times = T.times

using CairoMakie

n = Observable(1)

T_lims = (8.94, 9.06)
N_lims = (0, 4.5)
P_lims = (0.007, 0.02)

Tₙ = @lift interior(T[$n], :, 1, :)
Nₙ = @lift interior(N[$n], :, 1, :)
Pₙ = @lift interior(P[$n], :, 1, :)

fig = Figure(size = (1000, 520), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 770, yticks = [-400, -200, 0])
ax1 = Axis(fig[1, 1]; title = "Temperature (°C)", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "Nutrients concentration (mmol N / m³)",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, Tₙ, colorrange = T_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, Nₙ, colorrange = N_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [8.95, 9.0, 9.05])
Colorbar(fig[2, 2], hm2, ticks = [0, 2, 4])
Colorbar(fig[3, 2], hm3, ticks = [0.01, 0.02, 0.03])

rowgap!(fig.layout, 0)

record(fig, "buoyancy_front.gif", 1:length(times)) do i
    n[] = i
end