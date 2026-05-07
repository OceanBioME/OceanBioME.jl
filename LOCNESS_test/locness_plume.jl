import Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Random
using Printf
using NetCDF
using Statistics
using CUDA

using Oceananigans
using Oceananigans.Units: minute, minutes, hour, hours, day, days
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans: Forcing
using JLD2

using OceanBioME

# --- Constants & Data Setup ---
const total_moles_alk = 0.0e5 
const ρ₀ = 1026.0
const cp = 3991.0
const k_w = 2π / 60.0

# Interpolation helper (isbits compatible)
@inline function interpolate(t, y, T)
    if t <= T[1]
        return y[1]
    elseif t >= T[end]
        return y[end]
    else
        N = length(T)
        i = 2
        while i < N && T[i] < t
            i += 1
        end
        T₀, T₁ = T[i-1], T[i]
        y₀, y₁ = y[i-1], y[i]
        return y₀ + (y₁ - y₀) * (t - T₀) / (T₁ - T₀)
    end
end

# Profiling and initial state helpers
function LinearInterpolation(y, x)
    return z -> begin
        if z >= x[1] return y[1] end
        if z <= x[end] return y[end] end
        idx = findfirst(xi -> xi <= z, x)
        if isnothing(idx) return y[end] end
        if idx == 1 return y[1] end
        x1, x2 = x[idx-1], x[idx]
        y1, y2 = y[idx-1], y[idx]
        return y1 + (y2 - y1) * (z - x1) / (x2 - x1)
    end
end

function get_initial_profiles()
    # Use index 1 to match the start of the forcing (2023-09-01)
    idx = 1
    z_prof = -ncread("locness_thetao.nc", "depth")
    T_prof = ncread("locness_thetao.nc", "thetao")[1, 1, :, idx]
    S_prof = ncread("locness_so.nc", "so")[1, 1, :, idx]
    DIC_prof = ncread("locness_bgc.nc", "dissic")[1, 1, :, idx] # mol m⁻³
    Alk_prof = ncread("locness_bgc.nc", "talk")[1, 1, :, idx]   # mol m⁻³
    
    z_prof_uv = -ncread("locness_uo_vo.nc", "depth")
    U_prof = ncread("locness_uo_vo.nc", "uo")[1, 1, :, idx]
    V_prof = ncread("locness_uo_vo.nc", "vo")[1, 1, :, idx]

    # Handle NetCDF fill values and NaNs
    function clean_profile(p, default_val=0.0)
        p = [ (val > 1e10 || isnan(val)) ? NaN : val for val in p]
        if all(isnan.(p)) return fill(default_val, length(p)) end
        first_valid = p[findfirst(!isnan, p)]
        for i in 1:length(p)
            if isnan(p[i]) p[i] = first_valid end
        end
        return Float64.(p)
    end

    U_prof = clean_profile(U_prof, 0.0)
    V_prof = clean_profile(V_prof, 0.0)
    T_prof = clean_profile(T_prof, 20.0)
    S_prof = clean_profile(S_prof, 35.0)
    DIC_prof = clean_profile(DIC_prof, 2.0)
    Alk_prof = clean_profile(Alk_prof, 2.3)

    DIC_prof .*= 1000.0 # to mmol m⁻³
    Alk_prof .*= 1000.0 # to mmol m⁻³
    
    return z_prof, T_prof, S_prof, DIC_prof, Alk_prof, z_prof_uv, U_prof, V_prof
end

z_prof, T_prof, S_prof, DIC_prof, Alk_prof, z_prof_uv, U_prof, V_prof = get_initial_profiles()

T_int = LinearInterpolation(T_prof, z_prof)
S_int = LinearInterpolation(S_prof, z_prof)
DIC_int = LinearInterpolation(DIC_prof, z_prof)
Alk_int = LinearInterpolation(Alk_prof, z_prof)
U_int = LinearInterpolation(U_prof, z_prof_uv)
V_int = LinearInterpolation(V_prof, z_prof_uv)

# --- Grid ---
Nx = 128
Ny = 128
Nz = 40
Lx = 1000.0
Ly = 1000.0
Lz = 40.0

arch = CPU()
grid = RectilinearGrid(arch, size = (Nx, Ny, Nz),
                          x = (-Lx/2, Lx/2),
                          y = (-Ly/2, Ly/2),
                          z = (-Lz, 0),
                          topology = (Periodic, Periodic, Bounded))

# --- Forcing ---
function read_forcings(arch)
    # 1. Stokes drift
    vsdx_raw = ncread("locness_wav.nc", "VSDX")
    vsdy_raw = ncread("locness_wav.nc", "VSDY")
    sf_x = ncgetatt("locness_wav.nc", "VSDX", "scale_factor")
    off_x = ncgetatt("locness_wav.nc", "VSDX", "add_offset")
    sf_y = ncgetatt("locness_wav.nc", "VSDY", "scale_factor")
    off_y = ncgetatt("locness_wav.nc", "VSDY", "add_offset")
    
    vsdx = vsdx_raw .* sf_x .+ off_x
    vsdy = vsdy_raw .* sf_y .+ off_y

    time_wav_raw = ncread("locness_wav.nc", "time")
    uˢ₀_ts = Float64.(dropdims(mean(vsdx, dims=(1, 2)), dims=(1, 2)))
    vˢ₀_ts = Float64.(dropdims(mean(vsdy, dims=(1, 2)), dims=(1, 2)))
    t_wav = Float64.((time_wav_raw .- time_wav_raw[1]) .* 3600.0)
    
    # 2. Wind Stress
    τx_raw = ncread("locness_era5_instant.nc", "iews")
    τy_raw = ncread("locness_era5_instant.nc", "inss")
    time_era5_inst = ncread("locness_era5_instant.nc", "valid_time")
    τx_ts = Float64.(dropdims(mean(τx_raw, dims=(1, 2)), dims=(1, 2)))
    τy_ts = Float64.(dropdims(mean(τy_raw, dims=(1, 2)), dims=(1, 2)))
    t_era5 = Float64.(time_era5_inst .- time_era5_inst[1])
    
    # 3. Heat Flux
    ssr = ncread("locness_era5_accum.nc", "ssr")
    str = ncread("locness_era5_accum.nc", "str")
    sshf = ncread("locness_era5_accum.nc", "sshf")
    slhf = ncread("locness_era5_accum.nc", "slhf")
    time_era5_acc = ncread("locness_era5_accum.nc", "valid_time")
    q_net_ts = Float64.((dropdims(mean(ssr, dims=(1, 2)), dims=(1, 2)) .+ 
                         dropdims(mean(str, dims=(1, 2)), dims=(1, 2)) .+ 
                         dropdims(mean(sshf, dims=(1, 2)), dims=(1, 2)) .+ 
                         dropdims(mean(slhf, dims=(1, 2)), dims=(1, 2))) ./ 3600.0)
    t_q = Float64.(time_era5_acc .- time_era5_acc[1])
    
    return (τx = on_architecture(arch, τx_ts), 
            τy = on_architecture(arch, τy_ts), 
            t_era5 = on_architecture(arch, t_era5),
            q_net = on_architecture(arch, q_net_ts), 
            t_q = on_architecture(arch, t_q),
            uˢ₀ = on_architecture(arch, uˢ₀_ts), 
            vˢ₀ = on_architecture(arch, vˢ₀_ts), 
            t_wav = on_architecture(arch, t_wav))
end

forcings = read_forcings(arch)

# Oceananigans convention: FluxBoundaryCondition(Q) at top means -ν ∂z u = Q
# Physical wind stress τ = ρ ν ∂z u => ν ∂z u = τ/ρ => Q = -τ/ρ
@inline top_u_bc(x, y, t, p) = -interpolate(t, p.τx, p.t_era5) / ρ₀
@inline top_v_bc(x, y, t, p) = -interpolate(t, p.τy, p.t_era5) / ρ₀
# Similarly for heat: -κ ∂z T = Q. Heating (Q_net > 0) means Q < 0.
@inline top_T_bc(x, y, t, p) = -interpolate(t, p.q_net, p.t_q) / (ρ₀ * cp)

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(top_u_bc, parameters=(τx=forcings.τx, t_era5=forcings.t_era5)))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(top_v_bc, parameters=(τy=forcings.τy, t_era5=forcings.t_era5)))
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(top_T_bc, parameters=(q_net=forcings.q_net, t_q=forcings.t_q)))

@inline ∂z_uˢ(z, t, p) = 2 * k_w * interpolate(t, p.uˢ₀, p.t_wav) * exp(2 * k_w * z)
@inline ∂z_vˢ(z, t, p) = 2 * k_w * interpolate(t, p.vˢ₀, p.t_wav) * exp(2 * k_w * z)

stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ, ∂z_vˢ=∂z_vˢ, parameters=(uˢ₀=forcings.uˢ₀, vˢ₀=forcings.vˢ₀, t_wav=forcings.t_wav))

# --- Biogeochemistry ---
biogeochemistry = SimpleCaCO3Precipitation(; grid)

CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()
surface_CO2_flux = CenterField(grid)

# --- Virtual ship ---
t_deploy_start = 8hours
T_deploy = 75minutes
r_ship_max = 250.0
n_spiral_turns = 5.0
H_release = 5.0
fwhm_ship = 30.0

Q_alk_mmol_per_s = total_moles_alk * 1000.0 / T_deploy
spiral_pitch = r_ship_max / n_spiral_turns
σ_ship = fwhm_ship / (2 * sqrt(2 * log(2)))
alk_source_spatial_norm = 1 / (2π * σ_ship^2 * H_release)

@inline function ship_x(t, p)
    s = (t - p.t_start) / p.T_deploy
    θ = 2π * p.n_turns * s
    r = (p.spiral_pitch / (2π)) * θ
    return r * cos(θ)
end
@inline function ship_y(t, p)
    s = (t - p.t_start) / p.T_deploy
    θ = 2π * p.n_turns * s
    r = (p.spiral_pitch / (2π)) * θ
    return r * sin(θ)
end

@inline function alk_ship_forcing(x, y, z, t, p)
    (t < p.t_start || t > p.t_start + p.T_deploy) && return zero(x)
    (z <= -p.H_release || z >= 0) && return zero(x)
    xs = ship_x(t, p)
    ys = ship_y(t, p)
    dx = x - xs
    dy = y - ys
    gauss = exp(-(dx * dx + dy * dy) / (2 * p.σ^2))
    return p.Q_mmol_s * p.spatial_norm * gauss
end

alk_ship_params = (;
    t_start = t_deploy_start,
    T_deploy = T_deploy,
    n_turns = n_spiral_turns,
    spiral_pitch = spiral_pitch,
    σ = σ_ship,
    H_release = H_release,
    Q_mmol_s = Q_alk_mmol_per_s,
    spatial_norm = alk_source_spatial_norm,
)

alk_forcing = Forcing(alk_ship_forcing; parameters = alk_ship_params)

# --- Model ---
model = NonhydrostaticModel(grid;
                            advection = WENO(),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = SeawaterBuoyancy(),
                            biogeochemistry = biogeochemistry,
                            closure = AnisotropicMinimumDissipation(),
                            stokes_drift = stokes_drift,
                            boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, DIC = FieldBoundaryConditions(top = CO₂_flux)),
                            forcing = (Alk = alk_forcing,),
                            auxiliary_fields = (; surface_CO2_flux))

# --- Initial Conditions ---
Tᵢ(x, y, z) = T_int(z)
Sᵢ(x, y, z) = S_int(z)
DICᵢ(x, y, z) = DIC_int(z)
Alkᵢ(x, y, z) = Alk_int(z)

Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz)
uᵢ(x, y, z) = 1e-2 * Ξ(z) + U_int(z)
vᵢ(x, y, z) = 1e-2 * Ξ(z) + V_int(z)

set!(model, u=uᵢ, v=vᵢ, T=Tᵢ, S=Sᵢ, DIC=DICᵢ, Alk=Alkᵢ, CaCO₃=0.0)

# --- Simulation ---
stop_time = 36hours
simulation = Simulation(model, Δt=1.0, stop_time=stop_time)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=1minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

progress(sim) = @printf("Iteration: %d, time: %s, Δt: %s, max|w|: %.2e\n",
                        iteration(sim), prettytime(sim), prettytime(sim.Δt), maximum(abs, sim.model.velocities.w))

add_callback!(simulation, progress, IterationInterval(10))

function calculate_CO2_flux!(sim)
    clock = sim.model.clock
    grid = sim.model.grid
    Nz = grid.Nz
    for i in 1:grid.Nx
        for j in 1:grid.Ny
            CUDA.@allowscalar(sim.model.auxiliary_fields.surface_CO2_flux.data[i, j, Nz] = CO₂_flux.condition.func(i, j, grid, clock,
                (; DIC = sim.model.tracers.DIC,
                Alk = sim.model.tracers.Alk,
                T = sim.model.tracers.T,
                S = sim.model.tracers.S)))
            CUDA.@allowscalar(sim.model.auxiliary_fields.surface_CO2_flux.data[i, j, 1:Nz-1] .= sim.model.auxiliary_fields.surface_CO2_flux.data[i, j, Nz])
        end
    end
    return nothing
end
simulation.callbacks[:Calc_CO2_flux] = Callback(calculate_CO2_flux!, TimeInterval(4minutes))

# --- Output ---
filename = "locness"
simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers,
                                                filename = filename * "_tracers.jld2",
                                                schedule = TimeInterval(1hour),
                                                overwrite_existing = true)

simulation.output_writers[:slice_xz] =
    JLD2Writer(model, merge(model.velocities, model.tracers),
               filename = filename * "_xz.jld2",
               indices = (:, Int(grid.Ny/2), :),
               schedule = TimeInterval(4minutes),
               with_halos = false,
               overwrite_existing = true)

simulation.output_writers[:slice_xy] =
               JLD2Writer(model, merge(model.velocities, model.tracers),
                          filename = filename * "_xy.jld2",
                          indices = (:, :, grid.Nz),
                          schedule = TimeInterval(4minutes),
                          with_halos = false,
                          overwrite_existing = true)

simulation.output_writers[:CO2_flux] =
    JLD2Writer(model, (; surface_CO2_flux),
               filename = filename * "_CO2_flux.jld2",
               indices = (:, :, grid.Nz),
               schedule = TimeInterval(4minutes),
               with_halos = true,
               overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=TimeInterval(24hours), prefix="spinup_checkpoint")

run!(simulation)
