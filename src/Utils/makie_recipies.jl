using CairoMakie, OceanBioME, JLD2, RollingFunctions
using Oceananigans.Units

struct ResultsLOBSTER
    path
    t
    z
    NO₃
    NH₄
    P
    Z
    D
    DD
    Dᶜ
    DDᶜ
    DOM
    DIC
    ALK
    T
    S
    CO₂_flux
    sinking_flux
end

function ResultsLOBSTER(path, t, z, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, T, S, open_bottom::Bool = true)
    airsea, sinking_flux = zeros(length(t)), zeros(length(t))
    @show size(t)
    for (idx, time) in enumerate(t)
        airsea[idx] = Boundaries.airseaflux(0, 0, time, DIC[end, idx], ALK[end, idx], T(0, 0, 0, time), S(0, 0, 0, time), merge(Boundaries.defaults.airseaflux, (gas=:CO₂, )))
        sinking_flux[idx] = Dᶜ[1, idx]*LOBSTER.D_sinking + DDᶜ[1, idx]*LOBSTER.DD_sinking
    end
    return ResultsLOBSTER(path, t, z, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, T, S, airsea, sinking_flux)
end

function load_LOBSTER(path, T, S)
    file = jldopen(path)
    time_labels = keys(file["timeseries/NO₃"])[2:end]

    timeseries = (NO₃ = zeros(33, length(time_labels)), NH₄ = zeros(33, length(time_labels)), P = zeros(33, length(time_labels)), Z = zeros(33, length(time_labels)), D = zeros(33, length(time_labels)), DD = zeros(33, length(time_labels)), Dᶜ = zeros(33, length(time_labels)), DDᶜ = zeros(33, length(time_labels)), DOM = zeros(33, length(time_labels)), DIC = zeros(33, length(time_labels)), ALK = zeros(33, length(time_labels)), t=zeros(33, length(time_labels)))

    for (t, time) in enumerate(time_labels)
        for tracer in (LOBSTER.tracers..., :DIC, :ALK)
            timeseries[tracer][:, t] = file["timeseries/$tracer/$time"]
        end
        timeseries[:t][1, t] = file["timeseries/t/$time"]
    end

    z = file["grid/zᵃᵃᶜ"][4:end-3]
    return ResultsLOBSTER(path, timeseries.t[1, :], z, timeseries[(:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :DIC, :ALK)]..., T, S)
end

function Makie.plot!(res::ResultsLOBSTER)
    fig = Figure(resolution = (3600, 2000))
    fig[1, 1:8] = Label(fig, "LOBSTER results", tellwidth=false)

    ax_P = Axis(fig[2, 1], title = "Phytoplankton (mmol N/m³)")
    hm_P = heatmap!(ax_P, res.t/years, res.z, res.P')
    Colorbar(fig[2, 2], hm_P)

    ax_Z = Axis(fig[2, 3], title = "Zooplankton (mmol N/m³)")
    hm_Z = heatmap!(ax_Z, res.t/years, res.z, res.Z')
    Colorbar(fig[2, 4], hm_Z)

    ax_D = Axis(fig[2, 5], title = "Small detritus (mmol N/m³)")
    hm_D = heatmap!(ax_D, res.t/years, res.z, res.D')
    Colorbar(fig[2, 6], hm_D)

    ax_DD = Axis(fig[2, 7], title = "Large detritus (mmol N/m³)")
    hm_DD = heatmap!(ax_DD, res.t/years, res.z, res.DD')
    Colorbar(fig[2, 8], hm_DD)

    ax_DOM = Axis(fig[3, 1], title = "Disolved organic matter (mmol N/m³)")
    hm_DOM = heatmap!(ax_DOM, res.t/years, res.z, res.DOM')
    Colorbar(fig[3, 2], hm_DOM)

    ax_DIC = Axis(fig[3, 3], title = "Disolved inorganic matter (mmol N/m³)")
    hm_DIC = heatmap!(ax_DIC, res.t/years, res.z, res.DIC')
    Colorbar(fig[3, 4], hm_DIC)

    ax_ALK = Axis(fig[3, 5], title = "Alkalinity (mmol ⁻/m³)")
    hm_ALK = heatmap!(ax_ALK, res.t/years, res.z, res.ALK')
    Colorbar(fig[3, 6], hm_ALK)

    ax_flx = Axis(fig[3, 7:8], title = "Flux (mmol C/m²)")
    #airsea = lines!(ax_flx, res.t/years, res.CO₂_flux)
    airsea_mean = lines!(ax_flx, (res.t/years)[365:end], rollmean(res.CO₂_flux, 365), label="Air-sea CO₂ exchange")
    #sinking = lines!(ax_flx, res.t/years, res.sinking_flux)
    sinking_mean = lines!(ax_flx, (res.t/years)[365:end], rollmean(res.sinking_flux, 365), label="Sinking")
    axislegend(ax_flx; position = :rb)

    save("$(res.path[1:end-4])png", fig)

    @info "$(res.path[1:end-4])png"
    return fig
end