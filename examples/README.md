# Ocean Biogeochemial Modelling Environment Examples

## box.jl
In this example we will setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration. This demonstrates:
- How to setup OceanBioME's biogeochemical models as box models

## column.jl
In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstrates:
- How to setup OceanBioME's biogeochemical models
- How to setup light attenuation
- How to visualise results

This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first.

## kelp.jl
In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and active particles modelling the growth of sugar kelp. This demonstrates:
- How to setup OceanBioME's biogeochemical models
- How to setup light attenuation
- How to add biologically active particles which interact with the biodeochemical model
- How to include optional tracer sets (carbonate chemistry and oxygen)
- How to visualise results

This is forced by idealised mixing layer depth, surface photosynthetically available radiation (PAR) and temperature which are setup first.

## data_forced.jl
In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstrates:
- How to setup OceanBioME's biogeochemical models
- How to setup light attenuation
- How to load external forcing data
- How to run with optional tracer sets such as carbonate chemistry
- How to setup a non-uniform grid for better near surface resolution
- How to visualise results

The forcing data for this example is automatically downloaded before running, and details of its origin may be found in the NetCDF file. Forcing data is for a sub-polar area near Iceland with temperature and mixed layer depth data from [Copernicus](https://data.marine.copernicus.eu/viewer) and PAR data from [NASA](https://oceancolor.gsfc.nasa.gov/l3/).


## sediment.jl
In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model with a sediment model coupled providing the bottom boundary condition. This shows:
- How to setup OceanBioME's biogeochemical models
- How to setup light attenuation
- How to add a sediment (or other complicated) boundary model
- How to include optional tracer sets (carbonate chemistry and oxygen)
- How to visualise results
> WARNING: This sediment model is currently wrong and will not produce correct results
