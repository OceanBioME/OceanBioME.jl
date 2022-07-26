module Particles
#export of structs is frowned upon I think
using KernelAbstractions, Oceananigans, StructArrays
using Oceananigans.Architectures: device

@kernel function source_sink!(model, loc, source_properties, scalefactors, tracers, Δt)
    #this implimentation over doing it for single tracer at a time (to better match Oceananigans kernal function format) as identifying the particle mask takes time
    #don't like all the looping in this fucntion
    p = @index(Global)
    (fi, i::Int), (fj, j::Int), (fk, k::Int) = modf.(Oceananigans.Fields.fractional_indices(model.particles.properties.x[p], model.particles.properties.y[p], model.particles.properties.z[p], loc, model.grid))
   #possibly not the most efficient implimentation but benchmarked 10x faster than a weird vector solution        
    if (fi, fj, fk) != (0, 0, 0)
        d=zeros(2,2,2)
        vols=zeros(2,2,2)

        nodesi=zeros(Int, 2, 2, 2)
        nodesj=zeros(Int, 2, 2, 2)
        nodesk=zeros(Int, 2, 2, 2)
        for (a, di) in enumerate([fi, 1-fi]), (b, dj) in enumerate([fj, 1-fj]), (c, dk) in enumerate([fk, 1-fk])
            d[a,b,c] = sqrt(di^2+dj^2+dk^2)
            nodesi[a, b, c] = min(max(i+a, 1), model.grid.Nx)
            nodesj[a, b, c] = min(max(j+b, 1), model.grid.Ny)
            nodesk[a, b, c] = min(max(k+c, 1), model.grid.Nz)
        end

        normfactor=1/sum(1 ./ d)
        #have todo second loop to have the normfactor (I think)
        for (ind, i) in enumerate(nodesi)
            for (l, tracer) in enumerate(tracers)
                getproperty(model.tracers, tracer)[i, nodesj[ind], nodesk[ind]] += scalefactors[l] * Δt * normfactor * getproperty(model.particles.properties, source_properties[l])[p] / (d[ind] * Oceananigans.Operators.Vᶜᶜᶜ(i, nodesj[ind], nodesk[ind], model.grid))
            end
        end
    else
        for (l, tracer) in enumerate(tracers)
            getproperty(model.tracers, tracer)[i+1, j+1, k+1] += scalefactors[l] * Δt * getproperty(model.particles.properties, source_properties[l])[p] / (Oceananigans.Operators.Vᶜᶜᶜ(i, j, k, model.grid))
        end
    end
end

@kernel function property_update!(particles, equation::Function, arguments, params, update_properties, track_properties, Δt, t)
    p = @index(Global)

    @inbounds argument_values = [arg[p] for arg in arguments]

    results = equation(particles.properties.x, particles.properties.y, particles.properties.z, t, argument_values..., params, Δt)

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
        else
            tracer = getproperty(model.auxiliary_fields, tracked_field.tracer)
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
    property_update_event = property_update_kernal!(model.particles, model.particles.parameters.equation, function_arguments, function_parameters, model.particles.parameters.integrals, model.particles.parameters.diagnostics, Δt, model.clock.time)
    wait(property_update_event)

    #update source/sinks
    if length(particles.parameters.sink_fields)>0#probably not the most julian way todo this and would also be fixed if the next line was more generic
        LX, LY, LZ = location(getproperty(model.tracers, particles.parameters.sink_fields[1].tracer))#this will not be right if they are different between the differnet fields

        source_properties = getproperty.(particles.parameters.sink_fields, :property)
        source_tracers = getproperty.(particles.parameters.sink_fields, :tracer)
        source_scalefactors = getproperty.(particles.parameters.sink_fields, :scalefactor) .* model.particles.parameters.density

        sourcesink_kernal! = source_sink!(device(model.architecture), workgroup, worksize)
        sourcesink_event = sourcesink_kernal!(model, (LX(),LY(),LZ()), source_properties, source_scalefactors, source_tracers, Δt)
        wait(sourcesink_event)
    end
    particles.parameters.custom_dynamics(particles, model, Δt)
end

@inline no_dynamics(args...) = nothing

#Perhaps this should just be an overloading of the function name LagrangianParticles with different arguments
function setup(particles::StructArray, equation::Function, equation_arguments::NTuple{N, Symbol} where N, 
    equation_parameters::NamedTuple, integrals::NTuple{N, Symbol} where N, 
    diagnostics::NTuple{N, Symbol} where N, 
    tracked_fields::NTuple{N, NamedTuple{(:tracer, :property, :scalefactor), Tuple{Symbol, Symbol, Float64}}} where N, 
    sink_fields::NTuple{N, NamedTuple{(:tracer, :property, :scalefactor), Tuple{Symbol, Symbol, Float64}}} where N, 
    density::AbstractFloat, custom_dynamics=no_dynamics)
    
    for property in (equation_arguments..., integrals..., diagnostics...)
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
                                    diagnostics=diagnostics,
                                    density=density, 
                                    tracked_fields=tracked_fields,
                                    sink_fields=sink_fields,
                                    custom_dynamics=custom_dynamics
                                    ))
end
end#module
