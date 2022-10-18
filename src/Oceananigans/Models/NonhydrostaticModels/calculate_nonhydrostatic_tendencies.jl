using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device
using Oceananigans.Models.NonhydrostaticModels: calculate_Gu!, calculate_Gv!, calculate_Gw!, calculate_Gc!

using OceanBioME.Oceananigans.Fields: isacolumn
using OceanBioME.Oceananigans.NonhydrostaticModels: forced_auxiliary_fields

import Oceananigans.Models.NonhydrostaticModels: calculate_interior_tendency_contributions!

""" Store previous value of the source term and calculate current source term. """
function calculate_interior_tendency_contributions!(model)

    tendencies           = model.timestepper.Gⁿ
    arch                 = model.architecture
    grid                 = model.grid
    advection            = model.advection
    coriolis             = model.coriolis
    buoyancy             = model.buoyancy
    stokes_drift         = model.stokes_drift
    closure              = model.closure
    background_fields    = model.background_fields
    velocities           = model.velocities
    tracers              = model.tracers
    auxiliary_fields     = model.auxiliary_fields
    hydrostatic_pressure = model.pressures.pHY′
    diffusivities        = model.diffusivity_fields
    forcings             = model.forcing
    clock                = model.clock
    u_immersed_bc        = velocities.u.boundary_conditions.immersed
    v_immersed_bc        = velocities.v.boundary_conditions.immersed
    w_immersed_bc        = velocities.w.boundary_conditions.immersed

    workgroup, worksize = work_layout(grid, :xyz)

    calculate_Gu_kernel! = calculate_Gu!(device(arch), workgroup, worksize)
    calculate_Gv_kernel! = calculate_Gv!(device(arch), workgroup, worksize)
    calculate_Gw_kernel! = calculate_Gw!(device(arch), workgroup, worksize)
    calculate_Gc_kernel! = calculate_Gc!(device(arch), workgroup, worksize)

    barrier = Event(device(arch))


    Gu_event = calculate_Gu_kernel!(tendencies.u,
                                    grid,
                                    advection,
                                    coriolis,
                                    stokes_drift,
                                    closure,
                                    u_immersed_bc, 
                                    buoyancy,
                                    background_fields,
                                    velocities,
                                    tracers,
                                    auxiliary_fields,
                                    diffusivities,
                                    forcings,
                                    hydrostatic_pressure,
                                    clock,
                                    dependencies=barrier)

    Gv_event = calculate_Gv_kernel!(tendencies.v,
                                    grid,
                                    advection,
                                    coriolis,
                                    stokes_drift,
                                    closure,
                                    v_immersed_bc, 
                                    buoyancy,
                                    background_fields,
                                    velocities,
                                    tracers,
                                    auxiliary_fields,
                                    diffusivities,
                                    forcings,
                                    hydrostatic_pressure,
                                    clock,
                                    dependencies=barrier)

    Gw_event = calculate_Gw_kernel!(tendencies.w,
                                    grid,
                                    advection,
                                    coriolis,
                                    stokes_drift,
                                    closure,
                                    w_immersed_bc, 
                                    buoyancy,
                                    background_fields,
                                    velocities,
                                    tracers,
                                    auxiliary_fields,
                                    diffusivities,
                                    forcings,
                                    clock,
                                    dependencies=barrier)

    events = [Gu_event, Gv_event, Gw_event]

    for tracer_index in 1:length(tracers)
        @inbounds c_tendency = tendencies[tracer_index+3]
        @inbounds forcing = forcings[tracer_index+3]
        @inbounds c_immersed_bc = tracers[tracer_index].boundary_conditions.immersed

        Gc_event = calculate_Gc_kernel!(c_tendency,
                                        grid,
                                        Val(tracer_index),
                                        advection,
                                        closure,
                                        c_immersed_bc,
                                        buoyancy,
                                        background_fields,
                                        velocities,
                                        tracers,
                                        auxiliary_fields,
                                        diffusivities,
                                        forcing,
                                        clock,
                                        dependencies=barrier)

        push!(events, Gc_event)
    end

    for auxiliary_field in forced_auxiliary_fields(model)
        @inbounds a_tendency = tendencies[auxiliary_field]
        @inbounds forcing = forcings[auxiliary_field]

        workgroup, worksize = work_layout(grid, size(a_tendency))

        calculate_Ga_kernel! = calculate_Ga!(device(arch), workgroup, worksize)

        Ga_event = calculate_Ga_kernel!(a_tendency,
                                        grid,
                                        velocities,
                                        tracers,
                                        auxiliary_fields,
                                        forcing,
                                        clock,
                                        dependencies=barrier)

        push!(events, Ga_event)
    end

    wait(device(arch), MultiEvent(Tuple(events)))
    return nothing
end

@kernel function calculate_Ga!(Ga, args...)
    i, j, k = @index(Global, NTuple)
    @inbounds Ga[i, j, k] = auxiliary_tendency(i, j, k, args...)
end