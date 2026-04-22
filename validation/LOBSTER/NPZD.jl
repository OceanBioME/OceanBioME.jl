using OceanBioME, Oceananigans
using OceanBioME.Models.LOBSTERModel: Nutrient, Detritus

grid = RectilinearGrid(size = (1, 1, 10), extent = (1, 1, 10))

biogeochemistry = LOBSTER(; grid, biology = OceanBioME.Models.LOBSTERModel.PhytoZoo(temperature_coefficient = 1.88), nutrients = Nutrient(), detritus = Detritus(), carbonate_system = CarbonateSystem())

model = NonhydrostaticModel(grid; biogeochemistry, tracers = :T)

time_step!(model, 1)

# NPZ
biogeochemistry = LOBSTER(; grid, nutrients = Nutrient(), detritus = nothing, carbonate_system = CarbonateSystem())

model = NonhydrostaticModel(grid; biogeochemistry)

time_step!(model, 1)