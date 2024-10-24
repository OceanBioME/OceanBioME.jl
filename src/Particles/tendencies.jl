using Oceananigans: fields

function compute_particle_tendencies!(particles::BiogeochemicalParticles{N}, model) where N
    tendencies = particles.timestepper.tendencies

    field_names = required_particle_fields(particles)

    dev = device(architecture(model))

    workgroup = min(N, 256)
    worksize = N

    tendency_kernel! = _tendency_kernel!(dev, workgroup, worksize)

    for name in field_names
        tendency_kernel!(Val(name), particles, model.grid, fields(model), model.clock, tendencies)
    end

    return nothing
end

@kernel function _tendency_kernel!(val_field_name::Val{field_name}, particles, grid, fields, clock, tendencies) where field_name
    n = @index(Global)

    t = clock.time

    field_names = required_particle_fields(particles)
    Nf = length(field_names)

    field_values = ntuple(nf->particles.fields[nf][n], Val(Nf))

    tracer_values = extract_tracer_values(val_field_name, particles.field_interpolation, particles, grid, fields, n)

    tendencies[field_name][n] = particles.biogeochemistry(val_field_name, t, field_values..., tracer_values...)
end
