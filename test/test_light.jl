using OceanBioME, Test, Oceananigans

grid = RectilinearGrid(size=(1,1,2), extent=(1,1,2))

params=LOBSTER.defaults

@testset "Light attenuation" begin
    PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)

    model = NonhydrostaticModel(; grid, timestepper=:RungeKutta3, tracers=(:P, ), auxiliary_fields = (PAR=PAR, ))
    Pᵢ(x,y,z)=2.5+z
    set!(model, u=0, v=0, w=0, P=Pᵢ)
    sim = Simulation(model, Δt=0.1, stop_time=1)
    surface_PAR(t) = 1
    sim.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)));
    run!(sim)

    expected_PAR = [
        (exp(-params.k_r0*0.5-0.5*params.Χ_rp*(2*params.Rd_chl/params.r_pig)^params.e_r)*exp(-params.k_r0*1-1*params.Χ_rp*(1.5*params.Rd_chl/params.r_pig)^params.e_r)+exp(-params.k_b0*0.5-0.5*params.Χ_bp*(2*params.Rd_chl/params.r_pig)^params.e_b)*exp(-params.k_b0*1-1*params.Χ_bp*(1.5*params.Rd_chl/params.r_pig)^params.e_b))/2,
        (exp(-params.k_r0*0.5-0.5*params.Χ_rp*(2*params.Rd_chl/params.r_pig)^params.e_r)+exp(-params.k_b0*0.5-0.5*params.Χ_bp*(2*params.Rd_chl/params.r_pig)^params.e_b))/2
    ]

    results_PAR=convert(Array, model.auxiliary_fields.PAR)[1, 1, :]

    @test all(results_PAR .≈ expected_PAR)
end