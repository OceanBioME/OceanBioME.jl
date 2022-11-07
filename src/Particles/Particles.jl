module Particles
export ActiveLagrangianParticles
using KernelAbstractions, Oceananigans, StructArrays
using Oceananigans.Architectures: device, arch_array
using Oceananigans.Operators: Vᶜᶜᶜ
using Oceananigans.LagrangianParticleTracking: update_field_property!

include("tracer_tendencies.jl")

@kernel function solve_growth!(particles, equation::Function, arguments, params, prognostic, diagnostic, Δt, t)
    p = @index(Global)

    results = @inbounds equation(particles.properties.x[p], particles.properties.y[p], particles.properties.z[p], t, [arg[p] for arg in arguments]..., params, Δt)

    for property in prognostic
        @inbounds getproperty(particles.properties, property)[p] += getproperty(results, property) *Δt
    end

    for property in diagnostic
        @inbounds getproperty(particles.properties, property)[p] = getproperty(results, property)
    end
end

function particle_dynamics!(particles, model, Δt)
    num_particles = length(particles)
    workgroup = min(num_particles, 256)
    worksize = num_particles
    
    model_fields = fields(model)
    for (tracer_name, property_name) in pairs(particles.parameters.tracked_fields)
        tracer = model_fields[tracer_name]
        
        particle_property = getproperty(particles.properties, property_name)

        LX, LY, LZ = location(tracer)

        update_field_property_kernel! = update_field_property!(device(model.architecture), workgroup, worksize)
        source_event = update_field_property_kernel!(particle_property, particles.properties, model.grid, tracer, LX(), LY(), LZ())
        wait(source_event)
    end
    
    #update particle properties
    function_arguments = [getproperty(model.particles.properties, arg) for arg in model.particles.parameters.equation_arguments]
    function_parameters = model.particles.parameters.equation_parameters

    solve_growth_kernal! = solve_growth!(device(model.architecture), workgroup, worksize)
    solve_growth_event = solve_growth_kernal!(model.particles, model.particles.parameters.equation, function_arguments, function_parameters, model.particles.parameters.prognostic, model.particles.parameters.diagnostic, Δt, model.clock.time)
    wait(solve_growth_event)

    # execute user defined dynamics (e.g. floating)
    particles.parameters.custom_dynamics(particles, model, Δt)
end

function infinitesimal_particle_field_coupling!(model)
    num_particles = length(model.particles)
    workgroup = min(num_particles, 256)
    worksize = num_particles

    calculate_particle_tendency_kernel! =  calculate_particle_tendency!(device(model.architecture), workgroup, worksize)

    events = []
    for (tracer, property) in pairs(model.particles.parameters.coupled_fields)
        property_values = getproperty(model.particles.properties, property)
        tendency_field = model.timestepper.Gⁿ[tracer]

        calculate_particle_tendency_event = calculate_particle_tendency_kernel!(property_values, tendency_field, model.particles, model.grid)
        push!(events, calculate_particle_tendency_event)
    end
    
    wait(device(model.architecture), MultiEvent(Tuple(events)))
end

@inline no_dynamics(args...) = nothing

#Perhaps this should just be an overloading of the function name LagrangianParticles with different arguments
function ActiveLagrangianParticles(particles::StructArray;
    equation::Function = no_dynamics, 
    equation_arguments::NTuple{N, Symbol} where N = (), 
    equation_parameters::NamedTuple = NamedTuple(), 
    prognostic::NTuple{N, Symbol} where N = (), 
    diagnostic::NTuple{N, Symbol} where N = NamedTuple(), 
    tracked_fields::NamedTuple = NamedTuple(), 
    coupled_fields::NamedTuple = NamedTuple(),
    scalefactor::AbstractFloat = 1.0, 
    custom_dynamics=no_dynamics)
    
    required_fields = (equation_arguments..., prognostic..., diagnostic..., keys(tracked_fields)..., keys(coupled_fields)...)
    for property in required_fields
        if !(property in propertynames(particles))
            throw(ArgumentError("$property is a required particle property but $(eltype(particles)) has no property $property."))
        end
    end

    return LagrangianParticles(particles; 
                                        dynamics=particle_dynamics!, 
                                        parameters=(
                                            equation=equation,
                                            equation_arguments=equation_arguments,
                                            equation_parameters=equation_parameters,
                                            prognostic=prognostic,
                                            diagnostic=diagnostic,
                                            scalefactor=scalefactor, 
                                            tracked_fields=tracked_fields,
                                            coupled_fields=coupled_fields,
                                            custom_dynamics=custom_dynamics
                                            )
                                        )
end
end#module
