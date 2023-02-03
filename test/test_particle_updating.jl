using OceanBioME, Test, Oceananigans, StructArrays, CUDA
using OceanBioME.Particles: infinitesimal_particle_field_coupling!

struct CustomParticle{FT}
    #position
    x :: FT
    y :: FT
    z :: FT

    #properties
    A :: FT
    B :: FT

    #tracked fields
    C  :: FT
end

grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1), topology=(Periodic, Periodic, Periodic))

particlestruct = StructArray{CustomParticle}(([0.25], [0.25], [-0.25], [1.0], [0.0], [0.0]))

@testset "Passive particles" begin
    function particleupdate(x, y, z, t, A, B, params, Δt)
        return (A = -A, B = t)
    end

    particles = Particles.ActiveLagrangianParticles(particlestruct, 
        equation = particleupdate, 
        equation_arguments = (:A, :B),  
        prognostic = (:A, ),
        diagnostic = (:B, ),
        tracked_fields = (C = :C, ),
    )

    model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
    set!(model, u=0, v=0, w=0, C=1)
    sim = Simulation(model, Δt=0.1, stop_time=1)
    sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
    run!(sim)
    @testset "Tracer tracking" begin
        particle_c = convert(Array, model.particles.properties.C)
        model_c = convert(Array, model.tracers.C)
        @test size(particle_c) == tuple(1)
        @test all(particle_c .≈ model_c[1, 1, 1])
    end

    @testset "Property updating" begin
        particle_a = convert(Array, model.particles.properties.A)
        particle_b = convert(Array, model.particles.properties.B)

        @test particle_a[1] ≈ exp(-1) atol = 0.01
        @test all(particle_b .≈ 1.0)
    end
end

@testset "Active particles" begin
    function particleupdate(x, y, z, t, A, B, params, Δt)
        return (A=-0.1, B=t)
    end
    particlestruct=StructArray{CustomParticle}(([0.25], [0.25], [-0.25], [-0.1], [0.0], [0.0]))

    particles = Particles.ActiveLagrangianParticles(particlestruct, 
        equation = particleupdate, 
        equation_arguments = (:A, :B),  
        diagnostic = (:A, :B),
        tracked_fields = (C = :C, ),
        coupled_fields = (C = :A, )
    )

    model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
    set!(model, u=0, v=0, w=0, C=1)
    sim = Simulation(model, Δt=1, stop_time=1)
    sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
    run!(sim)

    CUDA.@allowscalar begin
        @test model.tracers.C[1, 1, 1] ≈ 0.9 atol = 0.05
    end

    @testset "Larger grid for point assignment" begin
        grid = RectilinearGrid(size=(2,1,1), extent=(2,1,1), topology=(Periodic, Periodic, Periodic))
        particlestruct=StructArray{CustomParticle}(([1.0], [0.25], [-0.25], [-0.1], [0.0], [0.0]))            

        particles = Particles.ActiveLagrangianParticles(particlestruct; 
            equation = particleupdate, 
            equation_arguments = (:A, :B),  
            diagnostic = (:A, :B), 
            tracked_fields = (C = :C, ),
            coupled_fields = (C = :A, )
        )

        model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
        set!(model, u=0, v=0, w=0, C=1)
        sim = Simulation(model, Δt=0.1, stop_time=1)
        sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
        run!(sim)

        CUDA.@allowscalar begin
            @test model.tracers.C[1, 1, 1] ≈ 0.95 atol = 0.05#twice the volume so should be half the conc. change
            @test model.tracers.C[2, 1, 1] == model.tracers.C[1, 1, 1] #evenly distributed so should spread equally
        end

        grid = RectilinearGrid(size=(4,4,4), extent=(4,4,4), topology=(Periodic, Periodic, Periodic))#have to use a bigger grid because (currently) there is an issue with BC enforecement in the particle trackign

        #on a grid point
        particlestruct=StructArray{CustomParticle}(([0.5], [0.5], [-1.5], [-0.1], [0.0], [0.0]))

        particles = Particles.ActiveLagrangianParticles(particlestruct; 
            equation = particleupdate, 
            equation_arguments = (:A, :B),  
            diagnostic = (:A, :B), 
            tracked_fields = (C = :C, ),
            coupled_fields = (C = :A, )
        )

        model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
        set!(model, u=0, v=0, w=0, C=1)
        sim = Simulation(model, Δt=1, stop_time=1)
        sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
        run!(sim)

        CUDA.@allowscalar begin
            @test model.tracers.C[1, 1, 3] ≈ 0.9
        end

        particlestruct=StructArray{CustomParticle}(([1.66], [1.55], [-1.74], [-0.1], [0.0], [0.0]))
        
        particles = Particles.ActiveLagrangianParticles(particlestruct; 
            equation = particleupdate, 
            equation_arguments = (:A, :B),  
            diagnostic = (:A, :B), 
            tracked_fields = (C = :C, ),
            coupled_fields = (C = :A, )
        )

        model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
        set!(model, u=0, v=0, w=0, C=1.0)
        sim = Simulation(model, Δt=1.0, stop_time=1)
        sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())

        run!(sim)

        CUDA.@allowscalar begin
            @test sum(model.tracers.C[1:4, 1:4, 1:4]) ≈ 4 ^ 3 - 0.1
        end

        #higher scalefactor
        particlestruct = StructArray{CustomParticle}(([1.66], [1.55], [-1.74], [-0.1], [0.0], [0.0]))

        function particleupdate(x, y, z, t, A, B, params, Δt)
            return (A=-0.1, B=t)
        end
        
        particles = Particles.ActiveLagrangianParticles(particlestruct; 
            equation = particleupdate, 
            equation_arguments = (:A, :B),  
            diagnostic = (:A, :B), 
            tracked_fields = (C = :C, ),
            coupled_fields = (C = :A, ),
            scalefactor = 10.0
        )

        model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, particles=particles, tracers=(:C, ))
        set!(model, u=0, v=0, w=0, C=1)
        sim = Simulation(model, Δt=1, stop_time=1)
        sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
                
        run!(sim)

        CUDA.@allowscalar  begin
            @test sum(model.tracers.C[1:4, 1:4, 1:4]) ≈ 4 ^ 3 - 1
        end
    end
end
