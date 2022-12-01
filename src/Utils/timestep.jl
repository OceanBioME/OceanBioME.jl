using Oceananigans.TurbulenceClosures: min_Δxyz, cell_diffusion_timescale, formulation, min_Δz
using Oceananigans.Units: day, minutes
using Oceananigans.BoundaryConditions: getbc
using Oceananigans: fields

#####
##### Find minimum diffusion timescale
#####

function minimum_timescale(κ::Function, grid, t::Number) 
    if !(grid.Nz==1 && grid.Ny ==1)
        Δz = isa(grid.Δzᵃᵃᶜ, Number) ? grid.Δzᵃᵃᶜ : grid.Δzᵃᵃᶜ[1:grid.Nz]
        return findmin([1/κ(0.5, 0.5, z, t) for z in grid.zᵃᵃᶜ[1:grid.Nz]].*Δz.^2)[1]
    else
        @warn "Multidimensional grid may take some time to find max functional diffusivity"
        diffs = zeros(grid.Nx, grid.Ny, grid.Nz)
        for (i, x) in enumerate(grid.xᶜᵃᵃ[1:grid.Nx]), (j, y) in enumerate(grid.yᵃᶜᵃ[1:grid.Ny]), (k, z) in enumerate(grid.zᵃᵃᶜ[1:grid.Nz])
            diffs[i, j, k] = min(grid.Δxᶜᵃᵃ[i], grid.Δyᵃᶜᵃ[j], grid.Δzᵃᵃᶜ[k])^2/κ(0.5, 0.5, z, t) 
        end
        return findmin(diffs)[1]
    end
end

minimum_timescale(κ::Function, grid, t::Vector)  = findmin([minimum_timescale(κ, grid, _t) for _t in t])[1]

function diffusion_timescale(model, grid; t=[1.0:365.0;].*day)
    if !isa(model.closure.κ[1], Function)
        return cell_diffusion_timescale(model.closure, nothing, grid)
    else
        min_κ = findmax([minimum_timescale(getproperty(model.closure.κ, tracer), grid, t) for tracer in keys(model.closure.κ)])[1]
        min_ν = minimum_timescale(model.closure.ν, grid, t)
        return min(min_κ, min_ν)
    end
end

#####
##### Find forcing error timescale
#####
function forcing_error_timescale(model)
    min_forcing = Inf
    
    for (tracer_idx, tracer) in enumerate(model.tracers)
        Gⁿ = @inbounds model.timestepper.Gⁿ[tracer_idx+3][1:model.grid.Nx, 1:model.grid.Ny, 1:model.grid.Nz]
        Gⁿ⁻¹ = @inbounds model.timestepper.G⁻[tracer_idx+3][1:model.grid.Nx, 1:model.grid.Ny, 1:model.grid.Nz]
        τⁿ⁺¹ = abs.(Gⁿ .- Gⁿ⁻¹)
        y = @inbounds tracer[1:model.grid.Nx, 1:model.grid.Ny, 1:model.grid.Nz]
        Δt = sqrt.(abs.(y./τⁿ⁺¹))
        Δt = Δt[isfinite.(Δt)]
        if length(Δt) > 0
            tracer_min = @inbounds findmin(Δt)[1]
            min_forcing =  min(min_forcing, tracer_min)
        end
    end
    return min_forcing
end

"""
    update_timestep!(sim, params=(w=200/day, c_diff = 0.45, c_adv = 0.45, relaxation=0.75))

Adapts the timestep (Δt) of `sim` to maintain stability.

Parameters
========
* `w`: maximum velocity for AdvectiveForcing
* `c_adv`: advective Courant number
* `c_diff`: diffusive Courant number
* `c_forcing`: forcing Courant number
* `Δt_max`: maximum timestep
* `relaxation`: exponent damping where Δtⁱ⁺¹ = Δtⁱ*(Cₘₐₓ/C)^relaxation
"""
function update_timestep!(sim, params=(w=200/day, c_diff = 0.45, c_adv = 0.45, relaxation=0.75))
    Δt_diff = :c_diff in keys(params) ? sim.Δt^(1-params.relaxation)*(params.c_diff*OceanBioME.diffusion_timescale(sim.model, sim.model.grid; t=sim.model.clock.time))^params.relaxation : Inf
    Δt_adv = :c_adv in keys(params) ? sim.Δt^(1-params.relaxation)*(params.c_adv/(abs(params.w)/OceanBioME.min_Δz(sim.model.grid)))^params.relaxation : Inf #replace with some way to check the actual max sinking velocity
    Δt_forcing = :c_forcing in keys(params) ? sim.Δt^(1-params.relaxation) * (params.c_forcing*forcing_error_timescale(sim.model))^params.relaxation : Inf
    Δt_max = :Δt_max in keys(params) ? params.Δt_max : Inf

    sim.Δt = findmin([Δt_diff, Δt_adv, Δt_forcing, Δt_max])[1]

    #Record the different timesteps for diagnostics if possible
    if :timesteps in keys(sim.model.auxiliary_fields) 
        sim.model.auxiliary_fields.timesteps[1] = [Δt_diff, Δt_adv, Δt_forcing]
    end
end
