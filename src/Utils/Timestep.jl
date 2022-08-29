using Oceananigans.TurbulenceClosures: min_Δxyz, cell_diffusion_timescale, formulation, min_Δz
using Oceananigans.Units: day, minutes
using Oceananigans.BoundaryConditions: getbc

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

function max_c_boundary(model, Δt)
    #only implimenting for top and bottom bcs for
    flux_ratios = []
    model_fields = merge(model.velocities, model.tracers)
    for (tracer_idx, tracer) in enumerate(model.tracers)
        for side in [:top, :bottom]
            boundary = getproperty(tracer.boundary_conditions, side)
            k = side==:top ? model.grid.Nz : 1
            if !isa(boundary.condition, Number) && !isa(boundary.condition, Nothing)
                for i=1:model.grid.Nx, j=1:model.grid.Ny
                    flux = Δt*getbc(boundary, i, j, model.grid, model.clock, model_fields)
                    value = tracer[i, j, k]
                    push!(flux_ratios, value!=0 ? flux/value : 0)
                end
            end
        end
    end
    return length(flux_ratios)>0 ? findmax(flux_ratios)[1] : 0
end

function update_timestep!(sim, params=(w=200/day, Δt_max=10minutes, c_diff = 0.6, c_adv = 0.8, relaxation=0.95))#, c_boundary=0.01))#off by default as computationally intensive and not limiting (as far as I can see)
    #possibly rewrite these as kernel functions if we're sticking with this kind of thing
    Δt_diff = :c_diff in keys(params) ? sim.Δt^(1-params.relaxation)*(params.c_diff*OceanBioME.diffusion_timescale(sim.model, sim.model.grid; t=sim.model.clock.time))^params.relaxation : Inf
    Δt_adv = :c_adv in keys(params) ? sim.Δt^(1-params.relaxation)*(params.c_adv/(abs(params.w)/OceanBioME.min_Δz(sim.model.grid)))^params.relaxation : Inf #replace with some way to check the actual max sinking velocity
    Δt_boundary = :c_boundary in keys(params) ? sim.Δt * (params.c_boundary/max_c_boundary(sim.model, sim.Δt))^params.relaxation : Inf
    Δt_max = :Δt_max in keys(params) ? params.Δt_max*(sim.Δt/params.Δt_max)^params.relaxation : Inf

    sim.Δt = min(Δt_diff, Δt_adv, Δt_boundary, Δt_max)
end