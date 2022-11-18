module Particles
export ActiveLagrangianParticles
using KernelAbstractions, Oceananigans, StructArrays
using Oceananigans.Architectures: device, arch_array
using Oceananigans.Operators: Vᶜᶜᶜ
using Oceananigans.LagrangianParticleTracking: update_field_property!

include("tracer_tendencies.jl")

# some how this is much faster for any length(args), but wrote because ntuple method is only fast for length(args)<11 as hard coded https://github.com/JuliaLang/julia/blob/master/base/ntuple.jl
@inline getargs(args, p) = length(args) < 11 ? ntuple(n->args[n][p], length(args)) : map(n -> arguments[n][p], 1:length(args))

@kernel function solve_growth!(particles, equation::Function, arguments, params, prognostic, diagnostic, Δt, t)
    p = @index(Global)
    @inbounds begin
        arg_values = getargs(arguments, p)
        results = equation(particles.properties.x[p], particles.properties.y[p], particles.properties.z[p], t, arg_values..., params, Δt)

        for property in prognostic
            getproperty(particles.properties, property)[p] += getproperty(results, property) *Δt
        end

        for property in diagnostic
            getproperty(particles.properties, property)[p] = getproperty(results, property)
        end
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

"""
    Particles.infinitesimal_particle_field_coupling!(model)

Applies tendencies specified for ActiveLagrangianParticles to fields in `model`.
Specifically with this coupling the tendencies are applied across the nearest nodes to the particle.

This should be used as a callback and the `TendencyCallsite` like:
```julia
sim.callbacks[:couple_particles] = Callback(Particles.infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
```
"""

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


"""
    Particles.ActiveLagrangianParticles(particles::StructArray;
                                                            equation::Function = no_dynamics, 
                                                            equation_arguments::NTuple{N, Symbol} where N = (), 
                                                            equation_parameters::NamedTuple = NamedTuple(), 
                                                            prognostic::NTuple{N, Symbol} where N = (), 
                                                            diagnostic::NTuple{N, Symbol} where N = NamedTuple(), 
                                                            tracked_fields::NamedTuple = NamedTuple(), 
                                                            coupled_fields::NamedTuple = NamedTuple(),
                                                            scalefactor::AbstractFloat = 1.0, 
                                                            custom_dynamics=no_dynamics)

Returns Oceananigans `LagrangianParticles` which will evolve and interact with biogeochemistry
Arguments
========
* `particles`: StructArray of the particles as per Oceananigans
Keyword arguments
=============
* `equation`: equation specifying how the particles properties and coupling evolves, should be of the form `equation(x, y, z, t, args..., params, Δt)` where args are specified by:
* `equation_arguments`: Tuple of Symbols specifying arguments for `equation`, arguments must be properties of `particles`
* `equation_parameters`: NamedTuple of parameters to pass to `equation`
* `prognostic`: properties for which `equation` returns dP/dt
* `diagnostic`: properties for which `equation` returns the new value (and therefore does not need to be integrarted)
* `tracked_fields`: NamedTuple specifying model fields to track where keys are the name of the field and values are names of the particle properties to store the fields in
* `coupled_fields`: NamedTuple specifying model fields coupled to the particles where keys are the field name and values are the property name to modify the fields *tendency* (i.e. the property should be a rate)
* `scalefactor`: multiplier to apply to field coupling, can be thought of as number of individuals that each particle represents
* `custom_dynamics`: custom function to apply to particles after `equation`, should be of the form `dynamics!(particles, model, Δt)`
```
"""

# Perhaps this should just be an overloading of the function name LagrangianParticles with different arguments
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
    
    required_fields = (equation_arguments..., prognostic..., diagnostic..., keys(tracked_fields)..., values(coupled_fields)...)
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
