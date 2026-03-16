using OceanBioME, Oceananigans
using OceanBioME.Models.LOBSTERModel: Nutrient, Detritus

grid = RectilinearGrid(size = (1, 1, 10), extent = (1, 1, 10))

biogeochemistry = LOBSTER(; grid, nutrients = Nutrient(), detritus = Detritus(), carbonate_system = CarbonateSystem(), oxygen = Oxygen())

model = NonhydrostaticModel(; grid, biogeochemistry)