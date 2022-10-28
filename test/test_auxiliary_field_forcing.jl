using Oceananigans, OceanBioME, Test, Statistics

grid = RectilinearGrid(size=(2, 2, 3), extent=(1, 1, 1))

B=Oceananigans.Fields.Field{Center, Center, Center}(grid)
C=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1)) 

for timestepper in (:QuasiAdamsBashforth2,  :RungeKutta3)
    @testset "Constant Forcing with $timestepper" begin
        forcing1(x, y, z, t) = 1.0

        model = NonhydrostaticModel(; grid, timestepper=timestepper, tracers=:A, auxiliary_fields=(; B, C), forcing = (B=forcing1, C=forcing1))
        set!(model, B=0.0, C=0.0)
        time_step!(model, 1)

        @test all(model.auxiliary_fields.B[1:2, 1:2, 1:3] .== 1.0)
        @test all(model.auxiliary_fields.C[1:2, 1:2, 1] .== 1.0)
    end

    @testset "Continuous tracer dependency with $timestepper" begin
        forcing2(x, y, z, t, A) = A
        Forcing2 = Forcing(forcing2, field_dependencies=:A)
        a₀(x, y, z) = x

        model = NonhydrostaticModel(; grid, timestepper=timestepper, tracers=:A, auxiliary_fields=(; B, C), forcing = (B=Forcing2, C=Forcing2))
        set!(model, A=a₀, B=0.0, C=0.0)

        time_step!(model, 1)
        @test all(model.auxiliary_fields.B[1:2, 1:2, 1:3] .== model.tracers.A[1:2, 1:2, 1:3])
        @test all(model.auxiliary_fields.C[1:2, 1:2, 1] .== model.tracers.A[1:2, 1:2, 1])
    end

    @testset "Continuous column integrated field dependency with $timestepper" begin
        #interpolations is modifed such that it is understood if a forcing has a dependency on a field that is Nx*Ny*1 
        #that it is a column integrated variable and as such returns the same value for every z
        forcing3(x, y, z, t, C) = C
        Forcing3 = Forcing(forcing3, field_dependencies=:C)
        c₀(x, y, z) = x

        model = NonhydrostaticModel(; grid, timestepper=timestepper, tracers=:A, auxiliary_fields=(; B, C), forcing = (A=Forcing3, B=Forcing3))
        set!(model, A=0.0, B=0.0, C=c₀)

        time_step!(model, 1)
        @test all(model.tracers.A[1:2, 1:2, 1:3] .== repeat(model.auxiliary_fields.C[1:2, 1:2, 1], 1, 1, 3))
        @test all(model.auxiliary_fields.B[1:2, 1:2, 1:3] .== repeat(model.auxiliary_fields.C[1:2, 1:2, 1], 1, 1, 3))
    end

    @testset "Discrete forcing of column function (realistic forcing of column integrated field) with $timestepper" begin
        function forcing4(i, j, k, grid, clock, model_fields)
            Ā = mean(model_fields.A, dims=3)
            B̄ = mean(model_fields.B, dims=3)
            return Ā[i, j, 1] + B̄[i, j, 1]
        end
        Forcing4 = Forcing(forcing4, discrete_form=true)

        a₀(x, y, z) = x
        b₀(x, y, z) = y

        model = NonhydrostaticModel(; grid, timestepper=timestepper, tracers=:A, auxiliary_fields=(; B, C), forcing = (C=Forcing4, ))
        set!(model, A=a₀, B=b₀, C=0.0)

        time_step!(model, 1)

        @test all(model.auxiliary_fields.C .== (mean(model.tracers.A, dims=3) + mean(model.auxiliary_fields.B, dims=3))[1:2, 1:2, 1])
    end
end