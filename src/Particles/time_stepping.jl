"""
    ForwardEuler

Step particle biogeochemistry with a `ForwardEuler` methods with `Δt` from
the physical model substep.
"""
struct ForwardEuler{T}
    tendencies :: T

    function ForwardEuler(field_names, number, arch)
        tendencies = NamedTuple{field_names}(ntuple(n->on_architecture(arch, zeros(number)), Val(length(field_names))))

        return new{typeof(tendencies)}(tendencies)
    end
end

# a different timestepper might want to loop over these etc.
function time_step_particle_fields!(timestepper::ForwardEuler, particles, model, Δt)
    compute_particle_tendencies!(particles, model)

    step_particle_biogeochemistry!(timestepper, particles, model, Δt)

    return nothing
end

function step_particle_biogeochemistry!(timestepper::ForwardEuler, particles::BiogeochemicalParticles{N}, model, Δt) where N
    tendencies = timestepper.tendencies

    field_names = required_particle_fields(particles)

    dev = device(architecture(model))

    workgroup = min(N, 256)
    worksize = N

    step_kernel! = _euler_step!(dev, workgroup, worksize)

    for name in field_names
        step_kernel!(particles.fields[name], tendencies[name], Δt)
    end

    return nothing
end

@kernel function _euler_step!(field, tendency, Δt)
    n = @index(Global)

    @inbounds field[n] += tendency[n] * Δt
end