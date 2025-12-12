"""
    compute_buffer_variable(::Val{:buffer_variable}, biogeochemistry, tracers...)

Function each model needs to implement to compute a `:buffer_variable` in the buffer.

`tracers` contains value of all tracer fields declared by `required_tracers` method for
the model.
"""
function compute_buffer_variable end

@inline function fill_all_buffer_variables!(particles::BiogeochemicalParticles{N}, model) where N

    dev = device(architecture(model))

    workgroup = min(N, 256)
    worksize = N
    kernel! = fill_single_buffer_variable!(dev, workgroup, worksize)

    for variable in buffer_variables(particles.biogeochemistry)
        kernel!(particles.buffer[variable], Val(variable), particles, model.grid, fields(model))
    end

    return nothing
end

@kernel function fill_single_buffer_variable!(buffer, variable_name, particles, grid, fields)
    n = @index(Global) 
    tracer_values = extract_tracer_values(variable_name, particles.field_interpolation, particles, grid, fields, n)
    @inbounds buffer[n] = compute_buffer_variable(variable_name, particles.biogeochemistry, tracer_values...)
    nothing
end
