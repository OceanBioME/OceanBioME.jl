using Oceananigans.TurbulenceClosures: min_Δxyz, cell_diffusion_timescale, formulation, min_Δz
using Oceananigans.Units: day, minutes

function minimum_timescale(κ::Function, grid, t::Number) 
    if !(grid.Nz==1 && grid.Ny ==1)
        return findmin([1/κ(0.5, 0.5, z, t) for z in grid.zᵃᵃᶜ[1:grid.Nz]].*grid.Δzᵃᵃᶜ[1:grid.Nz].^2)[1]
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

function update_timestep(sim, params=(w=200/day, Δt_max=10minutes, c_diff = 2.3, c_adv = 1))
    sim.Δt = @show min(c_diff*OceanBioME.diffusion_timescale(sim.model, grid; t=sim.model.clock.time), c_adv*OceanBioME.min_Δz(grid)/params.w, params.Δt_max)#replace with some way to check the actual max sinking velocity
end