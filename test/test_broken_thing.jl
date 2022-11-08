using Oceananigans
n=10
grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200))
model=NonhydrostaticModel(;grid=grid)
sim=Simulation(model, Δt=1, stop_iteration=1)
sim.callbacks[:test]=Callback(disp_model; callsite=TendencyCallsite())
run!(sim)