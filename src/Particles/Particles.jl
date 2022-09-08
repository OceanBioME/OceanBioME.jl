module Particles
using KernelAbstractions, Oceananigans, StructArrays
using Oceananigans.Architectures: device, arch_array
using Oceananigans.Operators: Vᶜᶜᶜ

function getmask(i, j, k, fi, fj, fk, grid)
    nodesi=arch_array(grid.architecture, zeros(Int, 2, 2, 2))
    nodesj=arch_array(grid.architecture, zeros(Int, 2, 2, 2))
    nodesk=arch_array(grid.architecture, zeros(Int, 2, 2, 2))
    
    d=arch_array(grid.architecture, zeros(2,2,2))

    for (a, di) in enumerate([fi, 1-fi]), (b, dj) in enumerate([fj, 1-fj]), (c, dk) in enumerate([fk, 1-fk])
        d[a,b,c] = sqrt(di^2+dj^2+dk^2)
        nodesi[a, b, c] = min(max(i+a, 1), grid.Nx)
        nodesj[a, b, c] = min(max(j+b, 1), grid.Ny)
        nodesk[a, b, c] = min(max(k+c, 1), grid.Nz)
    end

    normfactor=1/sum(1 ./ d)

    d ./= normfactor
    return nodesi, nodesj, nodesk, d
end

function apply_sinks!(model, particles, p, i, j, k, d, Δt)
    for sink in particles.parameters.sink_fields
        value =particles.parameters.density * sink.scalefactor * Δt * getproperty(model.particles.properties, sink.property)[p] / (d * Vᶜᶜᶜ(i, j, k, model.grid))
        if abs(value) > 0
            res = value + getproperty(model.tracers, sink.tracer)[i, j, k]
            if res < 0
                @warn "Particle tried to take tracer below zero, conserving tracer (may cause adverse effects on particle dynamics)"
                getproperty(model.tracers, sink.tracer)[i, j, k] = 0
                getproperty(particles.properties, sink.fallback)[p] += fallback(particles, p, res, sink.fallback_scalefactor)*Vᶜᶜᶜ(i, j, k, model.grid)
            else
                getproperty(model.tracers, sink.tracer)[i, j, k] += value
                getproperty(model.timestepper.Gⁿ, sink.tracer)[i, j, k] += value
            end
        end
    end
end

fallback(particles, p, diff, fallback_scalefactor::Number) = diff*fallback_scalefactor

fallback(particles, p, diff, fallback_scalefactor::NamedTuple{(:property, :constant), Tuple{Symbol, Float64}}) = diff*fallback_scalefactor.constant/(getproperty(particles.properties, fallback_scalefactor.property)[p])

@kernel function source_sink!(model, loc, Δt)
    p = @index(Global)
    (fi, i::Int), (fj, j::Int), (fk, k::Int) = modf.(Oceananigans.Fields.fractional_indices(model.particles.properties.x[p], model.particles.properties.y[p], model.particles.properties.z[p], loc, model.grid))

    if (fi, fj, fk) != (0, 0, 0)
        nodesi, nodesj, nodesk, d = getmask(i, j, k, fi, fj, fk, model.grid)
        for (ind, i) in enumerate(nodesi)
            apply_sinks!(model, model.particles, p, i, nodesj[ind], nodesk[ind], d[ind], Δt)
        end
    else
        apply_sinks!(model, model.particles, p, i+1, j+1, k+1, 1, Δt)
    end
end

@kernel function property_update!(particles, equation::Function, arguments, params, update_properties, track_properties, Δt, t)
    p = @index(Global)

    @inbounds argument_values = [arg[p] for arg in arguments]

    results = equation(particles.properties.x[p], particles.properties.y[p], particles.properties.z[p], t, argument_values..., params, Δt)

    for property in update_properties
        @inbounds getproperty(particles.properties, property)[p] += getproperty(results, property) *Δt
    end

    for property in track_properties
        @inbounds getproperty(particles.properties, property)[p] = getproperty(results, property)
    end
end

function dynamics!(particles, model, Δt)
    num_particles = length(particles)
    workgroup = min(num_particles, 256)
    worksize = num_particles
    
    #update tracked fields, won't need this when Oceananigans is fixed (https://github.com/CliMA/Oceananigans.jl/pull/2662)
    for tracked_field in particles.parameters.tracked_fields
        if tracked_field.tracer in keys(model.tracers)
            tracer = getproperty(model.tracers, tracked_field.tracer)
        elseif tracked_field.tracer in keys(model.auxiliary_fields)
            tracer = getproperty(model.auxiliary_fields, tracked_field.tracer)
        elseif tracked_field.tracer in (:u, :v, :w)
            tracer = getproperty(model.velocities, tracked_field.tracer)
        else
            throw(ArgumentError("$(tracked_field.tracer) is not a model field that can be tracked"))
        end
        
        particle_property = getproperty(particles.properties, tracked_field.property)

        LX, LY, LZ = location(tracer)

        update_field_property_kernel! = Oceananigans.LagrangianParticleTracking.update_field_property!(device(model.architecture), workgroup, worksize)
        source_event = update_field_property_kernel!(particle_property, particles.properties, model.grid, tracer, LX(), LY(), LZ())
        wait(source_event)
    end
    #update particle properties
    function_arguments = [getproperty.(model.particles.properties, arg) for arg in model.particles.parameters.equation_arguments]
    function_parameters = model.particles.parameters.equation_parameters

    property_update_kernal! = property_update!(device(model.architecture), workgroup, worksize)
    property_update_event = property_update_kernal!(model.particles, model.particles.parameters.equation, function_arguments, function_parameters, model.particles.parameters.integrals, model.particles.parameters.tracked, Δt, model.clock.time)
    wait(property_update_event)

    #update source/sinks
    if length(particles.parameters.sink_fields)>0#probably not the most julian way todo this and would also be fixed if the next line was more generic
        LX, LY, LZ = location(getproperty(model.tracers, particles.parameters.sink_fields[1].tracer))#this will not be right if they are different between the differnet fields

        sourcesink_kernal! = source_sink!(device(model.architecture), workgroup, worksize)
        sourcesink_event = sourcesink_kernal!(model, (LX(),LY(),LZ()), Δt)

        wait(sourcesink_event)
    end
    particles.parameters.custom_dynamics(particles, model, Δt)
end

@inline no_dynamics(args...) = nothing

#Perhaps this should just be an overloading of the function name LagrangianParticles with different arguments
function setup(particles::StructArray, equation::Function, equation_arguments::NTuple{N, Symbol} where N, 
    equation_parameters::NamedTuple, integrals::NTuple{N, Symbol} where N, 
    tracked::NTuple{N, Symbol} where N, 
    tracked_fields::NTuple{N, NamedTuple{(:tracer, :property, :scalefactor), Tuple{Symbol, Symbol, Float64}}} where N, 
    sink_fields, #::NTuple{N, NamedTuple{(:tracer, :property, :scalefactor, :fallback, :fallback_scalefactor), Tuple{Symbol, Symbol, Float64, Symbol, NamedTuple}}} where N - got too complicated
    density::AbstractFloat, custom_dynamics=no_dynamics)
    
    for property in (equation_arguments..., integrals..., tracked...)
        if !(property in propertynames(particles))
            throw(ArgumentError("$property is a required field but $(eltype(particles)) has no property $property."))
        end
    end

    for field in (tracked_fields..., sink_fields...)
        if !(field.property in propertynames(particles))
            throw(ArgumentError("$(field.property) is a required field for field tracking or source/sinking but $(eltype(particles)) has no property $property."))
        end
    end

    #add checks (e.g. check that all tracked fields properties are properties of the particles)
    return LagrangianParticles(particles; 
                                dynamics=dynamics!, parameters=(
                                    equation=equation,
                                    equation_arguments=equation_arguments,
                                    equation_parameters=equation_parameters,
                                    integrals=integrals,
                                    tracked=tracked,
                                    density=density, 
                                    tracked_fields=tracked_fields,
                                    sink_fields=sink_fields,
                                    custom_dynamics=custom_dynamics
                                    ))
end
end#module
