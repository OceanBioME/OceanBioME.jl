function update_tendencies!(bgc, particles::BiogeochemicalParticles{N}, model) where N
    tendencies = particles.timestepper.tendencies

    field_names = coupled_tracers(particles)

    dev = device(architecture(model))

    workgroup = min(N, 256)
    worksize = N

    tendency_kernel! = _update_tracer_tendencies!(dev, workgroup, worksize)

    for name in field_names
        tendency_kernel!(Val(name), particles, model.grid, fields(model), model.clock, model.timestepper.Gⁿ, particles.scalefactors)
    end

    return nothing
end

@kernel function _update_tracer_tendencies!(val_field_name::Val{field_name}, particles, grid, fields, clock, tendencies, scalefactors) where field_name
    n = @index(Global)

    t = clock.time

    field_names = required_particle_fields(particles)
    nf = length(field_names)

    field_values = ntuple(n->particles.fields[n], Val(nf))

    tracer_values = extract_tracer_values(particles.field_interpolation, particles, grid, fields, n)
    
    particle_tendency = particles.biogeochemistry(val_field_name, t, field_values..., tracer_values...) 

    scalefactor = @inbounds scalefactors[n]

    total_tendency = @show scalefactor * particle_tendency

    apply_tracer_tendency!(particles.field_interpolation, particles, grid, total_tendency, tendencies[field_name], n)
end