function update_tendencies!(bgc, particles::BiogeochemicalParticles{N}, model) where N
    field_names = coupled_tracers(particles)

    dev = device(architecture(model))

    workgroup = min(N, 256)
    worksize = N

    tendency_kernel! = _update_tracer_tendencies!(dev, workgroup, worksize)

    for field_name in field_names
        target_name = possible_tuple_value(field_names, field_name) # coupled_tracers can return a NamedTuple 

        tendency_kernel!(Val(field_name), Val(target_name), particles, model.grid, fields(model), model.clock, model.timestepper.Gⁿ, particles.scalefactors)
    end

    return nothing
end

possible_tuple_value(field_names, field_name) = field_names[field_name]
possible_tuple_value(::Tuple, field_name) = field_name

@kernel function _update_tracer_tendencies!(val_field_name::Val{field_name}, ::Val{target_name}, 
                                            particles, 
                                            grid, 
                                            fields, 
                                            clock, 
                                            tendencies, 
                                            scalefactors) where {field_name, target_name}
    n = @index(Global)

    t = clock.time

    field_names = required_particle_fields(particles)
    Nf = length(field_names)

    field_values = ntuple(nf->particles.fields[nf][n], Val(Nf))

    tracer_values = extract_tracer_values(particles.field_interpolation, particles, grid, fields, n)

    particle_tendency = particles.biogeochemistry(val_field_name, t, field_values..., tracer_values...) 

    scalefactor = @inbounds scalefactors[n]

    total_tendency = scalefactor * particle_tendency

    apply_tracer_tendency!(particles.field_interpolation, particles, grid, total_tendency, tendencies[target_name], n)
end