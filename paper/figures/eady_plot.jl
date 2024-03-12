using Oceananigans, GLMakie, JLD2

Nx = 512

# load all the fields

field_path = "eady.jld2"

field_all = FieldTimeSeries(field_path, "P")

times = field_all.times
grid = field_all.grid

P = [Array(interior(field_all[n], 1:Int(Nx / 4), :, :)) for n in eachindex(times)]

field_all = FieldTimeSeries(field_path, "NO₃")

N = [Array(interior(field_all[n], Int(Nx / 4 + 1):Int(2 * Nx / 4), :, :)) for n in eachindex(times)]

field_all = FieldTimeSeries(field_path, "NH₄")

for n in eachindex(times)
    N[n] .+= Array(interior(field_all[n], Int(Nx / 4 + 1):Int(2 * Nx / 4), :, :))
end

field_all = FieldTimeSeries(field_path, "sPOM")

OC = [Array(interior(field_all[n], Int(2 * Nx / 4 + 1):Int(3 * Nx / 4), :, :)) for n in eachindex(times)]

field_all = FieldTimeSeries(field_path, "bPOM")

for n in eachindex(times)
    OC[n] .+= Array(interior(field_all[n], Int(Nx / 4 + 1):Int(2 * Nx / 4), :, :))
end

field_all = FieldTimeSeries(field_path, "DOM")

for n in eachindex(times)
    OC[n] .+= Array(interior(field_all[n], Int(Nx / 4 + 1):Int(2 * Nx / 4), :, :))
end

field_all = FieldTimeSeries(field_path, "DIC")

DIC = [Array(interior(field_all[n], Int(3 * Nx / 4 + 1):Int(4 * Nx / 4), :, :)) for n in eachindex(times)]

# load particles

file = jldopen("eady_particles.jld2")

iterations = keys(file["timeseries/t"])

x, y, z, A, N_kelp, C = ntuple(n -> ones(36, length(iterations)) .* NaN, 6)

for (idx, it) in enumerate(iterations)
    particles = file["timeseries/particles/$it"]

    x[:, idx] = particles.x
    y[:, idx] = particles.y
    z[:, idx] = particles.z
    A[:, idx] = particles.A
    N_kelp[:, idx] = particles.N
    C[:, idx] = particles.C
end

close(file)

##### 
##### Animation
#####

# setup the observables

n = Observable(1)

N_plt   =   @lift N[$n]
P_plt   =   @lift P[$n]
OC_plt  =  @lift OC[$n]
DIC_plt = @lift DIC[$n]

x_plt = @lift x[:, $n]
y_plt = @lift y[:, $n]
z_plt = @lift z[:, $n]
A_plt = @lift A[:, $n]


fig = Figure(size = (1600, 1000))

ax = Axis3(fig[1:4, 1:4], aspect = (1, 1, 0.28),
           xticks = [0, 1000], yticks = [0, 1000], zticks = [-140, 0],
           xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
           xgridvisible = false, ygridvisible = false, zgridvisible = false,
           xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
           protrusions = (50, 30, 30, 30))

vm1 = contour!(ax, xc[1:Int(Nx / 4)], yc, zc, N_plt, levels = 50, colormap = Reverse(:bamako))
vm2 = contour!(ax, xc[Int(Nx / 4 + 1):Int(2 * Nx / 4)], yc, zc, P_plt, levels = 50, colormap = Reverse(:batlow))
vm3 = contour!(ax, xc[Int(2 * Nx / 4 + 1):Int(3 * Nx / 4)], yc, zc, OC_plt, levels = 50, colormap = :lajolla)
vm4 = contour!(ax, xc[Int(3 * Nx / 4 + 1):Int(4 * Nx / 4)], yc, zc, DIC_plt, levels = 50, colormap = Reverse(:devon))

sc = scatter!(ax, x_plt, y_plt, z_plt, color = A_plt, colormap=:grayC)

txt = text!(ax,
            [Point3f(xc[Int(1 + (i - 1) * Nx / 4)], yc[1], zc[1]) for i in 1:4],
            text = ["Nutrients", "Phytoplankton", "Organic carbon", "Inorganic carbon"],
            rotation = [0 for i in 1:4],
            align = (:left, :top),
            fontsize = 35,
            markerspace = :data)

supertitle = Label(fig[0, 1:4], "t = 0.0, n = 1")

record(fig, "all.mp4", eachindex(times)[2:end], framerate = 10) do i
    n[] = i
    supertitle.text = "t = $(prettytime(times[i])), n = $i"
    println("$i of $(length(times))")
end

n[] = 39
fig

#####
##### Static Figure
#####
n = 37

lims[1] = (min(minimum(N_plt[:, :, end]), minimum(N_plt[1, :, :]), minimum(N_plt[:, 1, :])), max(maximum(N_plt[:, :, end]), maximum(N_plt[1, :, :]), maximum(N_plt[:, 1, :])))


fig = Figure(size = (1600, 1000))

ax = Axis3(fig[1:4, 1:4], aspect = (1, 1, 0.28),
           xticks = [0, 1000], yticks = [0, 1000], zticks = [-140, 0],
           xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
           xgridvisible = false, ygridvisible = false, zgridvisible = false,
           xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
           protrusions = (50, 5, 5, 5))


N_plt   =   N[n]
P_plt   =   P[n]
OC_plt  =  OC[n]
DIC_plt = DIC[n]

x_plt = x[:, n]
y_plt = y[:, n]
z_plt = z[:, n]
A_plt = A[:, n]

vm1 = contour!(ax, xc[1:Int(Nx / 4)], yc, zc, N_plt, levels = 50, colormap = Reverse(:bamako), colorrange = lims[1])
vm2 = contour!(ax, xc[Int(Nx / 4 + 1):Int(2 * Nx / 4)], yc, zc, P_plt, levels = 50, colormap = Reverse(:batlow))
vm3 = contour!(ax, xc[Int(2 * Nx / 4 + 1):Int(3 * Nx / 4)], yc, zc, OC_plt, levels = 50, colormap = :lajolla)
vm4 = contour!(ax, xc[Int(3 * Nx / 4 + 1):Int(4 * Nx / 4)], yc, zc, DIC_plt, levels = 50, colormap = Reverse(:devon))

sc = scatter!(ax, x_plt, y_plt, z_plt, color = A_plt, colormap=:grayC)

txt = text!(ax,
            [Point3f(xc[Int(1 + (i - 1) * Nx / 4)], yc[1], zc[1]) for i in 1:4],
            text = ["| Nutrients", "| Phytoplankton", "| Organic carbon", "| Inorganic carbon"],
            rotation = [0 for i in 1:4],
            align = (:left, :top),
            fontsize = 33,
            markerspace = :data)

Colorbar(fig[5, 1], limits = (minimum(N_plt), maximum(N_plt)), colormap = Reverse(:bamako), vertical = false, label = "Nutrients (mmol N / m³)")
Colorbar(fig[5, 2], limits = (minimum(P_plt), maximum(P_plt)), colormap = Reverse(:batlow), vertical = false, label = "Phytoplankton (mmol N / m³)")
Colorbar(fig[5, 3], limits = (minimum(OC_plt), maximum(OC_plt)), colormap = :lajolla, vertical = false, label = "Organic Carbon (mmol C / m³)")
Colorbar(fig[5, 4], limits = (minimum(DIC_plt), maximum(DIC_plt)), colormap = Reverse(:devon), vertical = false, label = "Inorganic Carbon (mmol C / m³)")

save("eady.png", fig)
