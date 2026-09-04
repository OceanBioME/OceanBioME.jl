using Oceananigans.Grids: znode
using Oceananigans.Operators: Δzᵃᵃᶠ, Δzᵃᵃᶜ

abstract type AbstractSingleBandExponentialLightAttenuation{BA, IN, CA, SP} end # bands, interface, cell average, and surface light

const AbstractLight{BA, IN, CA, SP} = AbstractSingleBandExponentialLightAttenuation{BA, IN, CA, SP}

# single band
@kernel function integrate_light_attenuation!(la::AbstractLight{1, Nothing}, 
                                              PAR, grid, clock, Chl, surface_PAR)

    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)

    @inbounds PAR[i, j, grid.Nz] = PAR⁰ * surface_center_attenuation(i, j, grid, la, clock, Chl)

    @inbounds for k in grid.Nz-1:-1:1
        PAR[i, j, k] = PAR[i, j, k+1] * center_to_center_attenuation(i, j, k, grid, la, clock, Chl)
    end
end

# recording interfaces
@kernel function integrate_light_attenuation!(la::AbstractLight{1}, 
                                              PAR, PARᵢ, grid, clock, Chl, surface_PAR)

    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)
    
    @inbounds PARᵢ[i, j, grid.Nz+1] = PAR⁰

    @inbounds for k in grid.Nz:-1:1
        eᵏᵈᶻ = face_to_face_attenuation(i, j, k, grid, la, clock, Chl)

        PARᵢ[i, j, k] = PARᵢ[i, j, k+1] * eᵏᵈᶻ
        PAR[i, j, k] = -PARᵢ[i, j, k+1] * (1 - eᵏᵈᶻ)/log(eᵏᵈᶻ)
    end
end

@inline surface_center_attenuation(i, j, grid, la::AbstractLight{1}, clock, Chl) =
    exp(attenuation(i, j, grid.Nz, grid, la, clock, Chl) * znode(i, j, grid.Nz, grid, Center(), Center(), Center()))

@inline center_to_center_attenuation(i, j, k, grid, la::AbstractLight{1}, clock, Chl) =
    (exp(-attenuation(i, j, k+1, grid, la, clock, Chl) * Δzᵃᵃᶜ(i, j, k+1, grid)/2) *
     exp(-attenuation(i, j, k, grid, la, clock, Chl) * Δzᵃᵃᶜ(i, j, k, grid)/2))

@inline face_to_face_attenuation(i, j, k, grid, la::AbstractLight{1}, clock, Chl) =
    exp(-attenuation(i, j, k, grid, la, clock, Chl) * Δzᵃᵃᶜ(i, j, k, grid))

# multiple implicit bands
@kernel function integrate_light_attenuation!(la::AbstractLight{N, Nothing}, 
                                              PAR, grid, clock, Chl, surface_PAR) where N

    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)

    z = znode(i, j, grid.Nz, grid, Center(), Center(), Center())
    K = surface_center_attenuation(i, j, grid, la, clock, Chl, z)

    @inbounds PAR[i, j, grid.Nz] = PAR⁰ * total_attenuation(K, la)

    @inbounds for k in grid.Nz-1:-1:1
        Δz1 = Δzᵃᵃᶜ(i, j, k+1, grid)
        Δz2 = Δzᵃᵃᶜ(i, j, k, grid)
        K = center_to_center_attenuation(i, j, k, grid, la, clock, Chl, K, Δz1, Δz2)
        PAR[i, j, k] = PAR⁰ * total_attenuation(K, la)
    end
end

# recording interfaces
@kernel function integrate_light_attenuation!(la::AbstractLight{N}, 
                                              PAR, PARᵢ, grid, clock, Chl, surface_PAR) where N
    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)
    
    @inbounds PARᵢ[i, j, grid.Nz+1] = PAR⁰

    K = (one(grid), one(grid))

    @inbounds for k in grid.Nz:-1:1
        Δz = Δzᵃᵃᶜ(i, j, k, grid)
        K = face_to_face_attenuation(i, j, k, grid, la, clock, Chl, K, Δz)
        
        PARᵢ[i, j, k] = PAR⁰ * total_attenuation(K, la)

        eᵏᵈᶻ = PARᵢ[i, j, k] / PARᵢ[i, j, k+1]

        PAR[i, j, k] = -PARᵢ[i, j, k+1] * (1 - eᵏᵈᶻ)/log(eᵏᵈᶻ) # I'm not convinced this is correct and that we don't have to average with each component first
    end
end

@inline weight(::AbstractLight{N}, n) where N = 1/N

@inline @generated function total_attenuation(K, la::AbstractLight{N}) where N
    combined = Expr(:block)
    push!(combined.args, :(total = zero(K[1])))
    for n in 1:N
        push!(combined.args, :(total += @inbounds weight(la, $n) * K[$n]))
    end
    push!(combined.args, :(return total))

    return combined
end

@inline @generated function surface_center_attenuation(i, j, grid, la::AbstractLight{N}, clock, Chl, z) where N
    args = [:(exp(z * attenuation(i, j, grid.Nz, grid, la, clock, Chl, Val($n)))) for n in 1:N]

    return Expr(:tuple, args...)
end

@inline @generated function center_to_center_attenuation(i, j, k, grid, 
                                                         la::AbstractLight{N}, 
                                                         clock, Chl, 
                                                         cumulative_attenuation, 
                                                         Δz1, Δz2) where N
    args = [:(@inbounds cumulative_attenuation[$n] * 
                        (exp(-Δz1/2 * attenuation(i, j, k+1, grid, la, clock, Chl, Val($n))) *
                         exp(-Δz2/2 * attenuation(i, j, k, grid, la, clock, Chl, Val($n))))) for n in 1:N]

    return Expr(:tuple, args...)
end

@inline @generated function face_to_face_attenuation(i, j, k, grid, la::AbstractLight{N}, clock, Chl, cumulative_attenuation, Δz) where N
    total_args = [:(@inbounds cumulative_attenuation[$n] * exp(-Δz * attenuation(i, j, k, grid, la, clock, Chl, Val($n)))) for n in 1:N]
    return Expr(:tuple, total_args...)
end

function update_biogeochemical_state!(model, PAR::AbstractLight{<:Any, Nothing})
    arch = architecture(model.grid)

    launch!(arch, model.grid, :xy, integrate_light_attenuation!, 
            PAR, PAR.field, model.grid, model.clock, 
            chlorophyll(model.biogeochemistry, model), 
            PAR.surface_PAR)

    return nothing
end

function update_biogeochemical_state!(model, PAR::AbstractLight)
    arch = architecture(model.grid)

    launch!(arch, model.grid, :xy, integrate_light_attenuation!, 
            PAR, PAR.field, PAR.interface_field, model.grid, 
            model.clock, chlorophyll(model.biogeochemistry, model), 
            PAR.surface_PAR)

    return nothing
end
