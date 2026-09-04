using Oceananigans.Operators: Δzᵃᵃᶜ

abstract type AbstractSingleBandExponentialLightAttenuation{BA, IN, CA, SP} end # bands, interface, cell average, and surface light

const AbstractLight{BA, IN, CA, SP} = AbstractSingleBandExponentialLightAttenuation{BA, IN, CA, SP}

# single band
@kernel function integrate_light_attenuation!(la::AbstractLight{1, Nothing}, 
                                              PAR, grid, clock, Chl, surface_PAR)

    i, j = @index(Global, NTuple)

    PARᵢ = getbc(surface_PAR, i, j, grid, clock, Chl)

    @inbounds for k in grid.Nz:-1:1
        eᵏᵈᶻ = face_to_face_attenuation(i, j, k, grid, la, clock, Chl)
        PAR[i, j, k] = - PARᵢ * (1 - eᵏᵈᶻ)/log(eᵏᵈᶻ)
        PARᵢ *= eᵏᵈᶻ
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

@inline face_to_face_attenuation(i, j, k, grid, la::AbstractLight{1}, clock, Chl) =
    exp(-attenuation(i, j, k, grid, la, clock, Chl) * Δzᵃᵃᶜ(i, j, k, grid))

# multiple implicit bands
@kernel function integrate_light_attenuation!(la::AbstractLight{N, Nothing}, 
                                              PAR, grid, clock, Chl, surface_PAR) where N

    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)
    K = initial_attenuation(la, eltype(grid))

    @inbounds for k in grid.Nz:-1:1
        Δz = Δzᵃᵃᶜ(i, j, k, grid)
        K_next = face_to_face_attenuation(i, j, k, grid, la, clock, Chl, K, Δz)

        PAR[i, j, k] = - PAR⁰ * total_cell_average(K, K_next, la)

        K = K_next
    end
end

# recording interfaces
@kernel function integrate_light_attenuation!(la::AbstractLight{N}, 
                                              PAR, PARᵢ, grid, clock, Chl, surface_PAR) where N
    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, Chl)
    
    @inbounds PARᵢ[i, j, grid.Nz+1] = PAR⁰

    K = initial_attenuation(la, eltype(grid))

    @inbounds for k in grid.Nz:-1:1
        Δz = Δzᵃᵃᶜ(i, j, k, grid)
        K_next = face_to_face_attenuation(i, j, k, grid, la, clock, Chl, K, Δz)

        PARᵢ[i, j, k] = PAR⁰ * total_attenuation(K_next, la)

        PAR[i, j, k] = - PAR⁰ * total_cell_average(K, K_next, la)

        K = K_next
    end
end

@inline initial_attenuation(::AbstractLight{N}, FT = Float64) where N = ntuple(_ -> one(FT), N)
@inline weight(::AbstractLight{N}, val_n) where N = 1/N

@inline @generated function total_attenuation(K, la::AbstractLight{N}) where N
    combined = Expr(:block)
    push!(combined.args, :(total = zero(K[1])))
    for n in 1:N
        push!(combined.args, :(total += @inbounds weight(la, Val($n)) * K[$n]))
    end
    push!(combined.args, :(return total))

    return combined
end

@inline @generated function face_to_face_attenuation(i, j, k, grid, la::AbstractLight{N}, clock, Chl, cumulative_attenuation, Δz) where N
    total_args = [:(@inbounds cumulative_attenuation[$n] * exp(-Δz * attenuation(i, j, k, grid, la, clock, Chl, Val($n)))) for n in 1:N]
    return Expr(:tuple, total_args...)
end

@inline @generated function total_cell_average(K, K_next, la::AbstractLight{N}) where N
    combined = Expr(:block)
    push!(combined.args, :(total = zero(K[1])))
    for n in 1:N
        push!(combined.args, :(total += @inbounds weight(la, Val($n)) * K[$n] * (1 - K_next[$n]/K[$n]) / log(K_next[$n]/K[$n])))
    end
    push!(combined.args, :(return total))

    return combined
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
