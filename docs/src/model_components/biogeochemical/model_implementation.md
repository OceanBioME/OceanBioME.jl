# [Implementing a new models](@id model_implementation)

Here we will describe how OceanBioME defines biogeochemical (BGC) models, how this varies from Oceananigans, and how to impliemnt your own model.

## Model structure
OceanBioME BGC models are `struct`s of type `ContinuousFormBiogeochemistry`, which is of abstract type `AbstractContinuousFormBiogeochemistry` from Oceananigans. In Oceananigans this describes BGC models which are continuous (depend continuously on ``x``, ``y``, and ``z``) rather than discrete (depending on ``i``, ``j``, ``k``). This simplifies the implementation of BGC models as forcing and sinking are very simple to define and add to tracers, and then Oceananigans handles the rest.

OceanBioME's `ContinuousFormBiogeochemistry` adds a layer on top of this which makes adding [light attenuation models](@ref light), [sediment](@ref sediment), and [biologically active particles](@ref individuals). This is because `ContinuousFormBiogeochemistry` has parameters in which the types of these components are stored. This means that the model components will automatically be integrated with the BGC model without having to add new methods various Oceananigans functions as you otherwise will have. 

## Implementing a model

The nature of multiple dispatch in Julia means that we define new BGC models as new types. You can then define [methods](https://docs.julialang.org/en/v1/manual/methods/) to this type which are used by OceanBioME and Oceananigans to integrate the model.

### The basics

For this example we are going to implement the simple Nutrient-Phytoplankton model similar to that used in [Chen2015](@citep), although we neglect the nutrient in/outflow terms since they may be added as [boundary conditions](https://clima.github.io/OceananigansDocumentation/stable/model_setup/boundary_conditions/), and modified to conserve nitrogen.

The first step we need is to import the abstract type from OceanBioME, some units from Oceananigans (for ease of parameter definition), and [`import`](https://stackoverflow.com/questions/27086159/what-is-the-difference-between-using-and-import-in-julia-when-building-a-mod) some functions from Oceananigans which we will add methods to:

```@example implementing
using OceanBioME: ContinuousFormBiogeochemistry
using Oceananigans.Units

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!,
                                     biogeochemical_drift_velocity,
				                     biogeochemical_auxiliary_fields
```

We then define our `struct` with the model parameters, as well as slots for the particles, light attenuation, and sediment models:

```@example implementing
@kwdef struct NutrientPhytoplankton{FT, LA, S, W, P} <:ContinuousFormBiogeochemistry{LA, S, P}
            base_growth_rate :: FT = 1.27 / day              # 1 / seconds
    nutrient_half_saturation :: FT = 0.025 * 1000 / 14       # mmol N / m³
       light_half_saturation :: FT = 300.0                   # micro einstein / m² / s
        temperature_exponent :: FT = 0.24                    # 1
     temperature_coefficient :: FT = 1.57                    # 1
         optimal_temperature :: FT = 28.0                    # °C
              mortality_rate :: FT = 0.15 / day              # 1 / seconds
     crowding_mortality_rate :: FT = 0.004 / day / 1000 * 14 # 1 / seconds / mmol N / m³

     light_attenuation_model :: LA = nothing
              sediment_model :: S  = nothing
            sinking_velocity :: W  = 2 / day
                   particles :: P  = nothing
end
```

In the above we used `@kwdef` to set default values for the models so that we don't have to set these up each time we use the model. We have also included a `sinking_velocity` field to later demonstrate how we can get tracers to sink. We also need to define some functions so that OceanBioME and Oceananigans know what tracers and auxiliary fields (e.g. light intensity) we need setting up:

```@example implementing
required_biogeochemical_tracers(::NutrientPhytoplankton) = (:N, :P, :T)

required_biogeochemical_auxiliary_fields(::NutrientPhytoplankton) = (:PAR, )

biogeochemical_auxiliary_fields(bgc::NutrientPhytoplankton) = biogeochemical_auxiliary_fields(bgc.light_attenuation_model)
```

Next we need to define the methods that specify how the nutrient, ``N``, and phytoplankton ``P`` evolve. We want the phytoplankton to evolve at the rate given by:

```math
\frac{\partial P}{\partial t} = \mu g(T) f(N) h(PAR) P - mP - bP^2,
```

where ``\mu`` is `base_growth_rate`, ``m`` is `mortality_rate`, ``b`` is `crowding`, and ``g``, ``f``, and ``h`` are given by:
```math
\begin{align}
g(T) &= c_1\exp\left(-c_2|T - T_{opt}|\right),\\
f(N) &= \frac{N}{k_N + N},\\
h(PAR) &= \frac{PAR}{k_P + PAR},
\end{align}
```
where ``c_1`` and ``c_2`` are `temperature_coefficient` and `temperature_exponent`, ``T_{opt}`` is `optimal_temperature`, ``k_N`` is `nutrient_half_saturation`, and ``k_P`` is `light_half_saturation`. Since this is a simple two variable model where the total concentration is conserved, ``\frac{\partial N}{\partial t} = - \frac{\partial P}{\partial t}``.

We turn this into a method of our new type by writing:

```@example implementing
@inline function (bgc::NutrientPhytoplankton)(::Val{:P}, x, y, z, t, N, P, T, PAR)
    μ = bgc.base_growth_rate
    m = bgc.mortality_rate
    b = bgc.crowding_mortality_rate

    growth = μ * g(bgc, T) * f(bgc, N) * h(bgc, PAR) * P

    death = m * P + b * P ^ 2

    return growth - death
end

@inline function g(bgc, T)
    c₁ = bgc.temperature_coefficient
    c₂ = bgc.temperature_exponent
    Tₒ = bgc.optimal_temperature

    return c₁ * exp(-c₂ * abs(T - Tₒ))
end

@inline function f(bgc, N)
    kₙ = bgc.nutrient_half_saturation

    return N / (N + kₙ)
end

@inline function h(bgc, PAR)
    kₚ = bgc.light_half_saturation

    return PAR / (PAR + kₚ)
end
```

The first parameter `::Val{:P}` is a special [value type](http://www.jlhub.com/julia/manual/en/function/Val) that allows this function to be dispatched when it is given the value `Val(:P)`. This is how Oceananigans tells the model which forcing function it wants. At the start of this function we unpack some parameters from the model, then calculate each term, and return the total change. For the nutrient we will define the evolution by:

```@example implementing
@inline (bgc::NutrientPhytoplankton)(::Val{:N}, args...) = -bgc(Val(:P), args...)
```

Finally, we need to define some functions to make sure that the sub-models (light attenuation) get updated in the right place

```@example implementing
using OceanBioME: BoxModel
import OceanBioME.BoxModels: update_boxmodel_state!

function update_biogeochemical_state!(bgc::NutrientPhytoplankton, model)
    update_PAR!(model, bgc.light_attenuation_model)
end

function update_boxmodel_state!(model::BoxModel{<:NutrientPhytoplankton, <:Any, <:Any, <:Any, <:Any, <:Any})
    getproperty(model.values, :PAR) .= model.forcing.PAR(model.clock.time)
    getproperty(model.values, :T) .= model.forcing.T(model.clock.time)
end
```

Now we can run an example similar to the [LOBSTER box model example](@ref box_example):

```@example implementing
using OceanBioME, Oceananigans.Units

const year = years = 365days

@inline PAR⁰(t) = 500 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # specify the nominal depth of the box for the PAR profile
@inline PAR(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay 

@inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

biogeochemistry = NutrientPhytoplankton() 

model = BoxModel(; biogeochemistry, forcing = (; PAR, T = temp))
model.Δt = 5minutes
model.stop_time = 5years

set!(model, N = 15, P = 15)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(model, save_interval = 100, save = SaveBoxModel("box_np.jld2"))
```

<details>
<summary>We can then visualise this:</summary>

```@example implementing
using JLD2

vars = (:N, :P, :T, :PAR)
file = jldopen("box_np.jld2")
times = parse.(Float64, keys(file["values"]))

timeseries = NamedTuple{vars}(ntuple(t -> zeros(length(times)), length(vars)))

for (idx, time) in enumerate(times)
    values = file["values/$time"]
    for tracer in vars
        getproperty(timeseries, tracer)[idx] = values[tracer]
    end
end

close(file)

# ## And plot
using CairoMakie

fig = Figure(resolution = (1200, 480), fontsize = 20)

axN= Axis(fig[1, 1], ylabel = "Nutrient \n(mmol N / m³)")
lines!(axN, times / year, timeseries.N, linewidth = 3)

axP = Axis(fig[1, 2], ylabel = "Phytoplankton \n(mmol N / m³)")
lines!(axP, times / year, timeseries.P, linewidth = 3)

axPAR= Axis(fig[2, 1], ylabel = "PAR (einstein / m² / s)", xlabel = "Time (years)")
lines!(axPAR, times / year, timeseries.PAR, linewidth = 3)

axT = Axis(fig[2, 2], ylabel = "Temperature (°C)", xlabel = "Time (years)")
lines!(axT, times / year, timeseries.T, linewidth = 3)

save("box_np.png", fig)
```
</details>

![buoyancy_front_NP](buoyancy_front_NP.mp4)

So now we know it works.

### Phytoplankton sinking

Now that we have a fully working BGC model we might want to add some more features. Another aspect that is easy to add is negative buoyancy. To-do this all we do is add a method to the Oceananigans function `biogeochemical_drift_velocity`. In this case we are going to only consider the phytoplankton to sink, so when we define the velocity we just need to set `sinking_velocity` to the speed like `NutrientPhytoplankton(; light_attenuation_model, sinking_velocity = 2/day)`. Then we add some more Oceananigans functions, and add the method:

```@example implementing
using Oceananigans.Fields: ZeroField, ConstantField

biogeochemical_drift_velocity(bgc::NutrientPhytoplankton, ::Val{:P}) = 
    (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocity)
```

### Sediment model coupling

Another aspect that OceanBioME allows easy coupling with is sediment models. Doing this varies between sediment models, but for the most generic and simplest all we need to do is add methods to two functions:

```@example implementing
using OceanBioME.Boundaries.Sediments: sinking_flux

import OceanBioME.Boundaries.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver

@inline nitrogen_flux(grid, advection, bgc::NutrientPhytoplankton, tracers, i, j) =
    sinking_flux(i, j, grid, advection, Val(:P), bgc, tracers)
                 
@inline carbon_flux(bgc::NutrientPhytoplankton, tracers, i, j) = nitrogen_flux(bgc, tracers, i, j) * 6.56

@inline remineralisation_receiver(::NutrientPhytoplankton) = :N
```

### Putting it together

Now that we have added these elements we can put it together into another simple example:
```@example implementing
using Oceananigans, OceanBioME
using OceanBioME.Sediments: InstantRemineralisation

# define some simple forcing

@inline surface_PAR(x, y, t) = 200 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 28

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - 50) / 10)) / 2 + 1e-4

# define the grid

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (32, ), extent = (100, ))

# setup the biogeochemical model

light_attenuation_model = TwoBandPhotosyntheticallyActiveRadiation(; grid, surface_PAR)

sediment_model = InstantRemineralisation(; grid, redfield = 6.56 * 14)

sinking_velocity = ZFaceField(grid)

w_sink(x, y, z) = 2 / day * tanh(z / 5)

set!(sinking_velocity, w_sink)

biogeochemistry = NutrientPhytoplankton(; light_attenuation_model, sinking_velocity, sediment_model) 

# put the model together

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry,
                              forcing = (T = Relaxation(rate = 1/day, target=temp), ))

set!(model, P = 0.01, N = 15, T = 28)

# run

simulation = Simulation(model, Δt = 9minutes, stop_time = 1years)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "column_np.jld2",
                                                       schedule = TimeInterval(1day),
                                                       overwrite_existing = true)

simulation.output_writers[:sediment] = JLD2OutputWriter(model, model.biogeochemistry.sediment_model.fields,
                                                        indices = (:, :, 1),
                                                        filename = "column_np_sediment.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:N, :P))
simulation.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())

run!(simulation)
```

<details>
<summary>We can then visualise this:</summary>

```@example implementing
N = FieldTimeSeries("column_np.jld2", "N")
P = FieldTimeSeries("column_np.jld2", "P")

sed = FieldTimeSeries("column_np_sediment.jld2", "N_storage")

fig = Figure()

axN = Axis(fig[1, 1], ylabel = "z (m)")
axP = Axis(fig[2, 1], ylabel = "z (m)")
axSed = Axis(fig[3, 1:2], ylabel = "Sediment (mmol N / m²)", xlabel = "Time (years)")

_, _, zc = nodes(grid, Center(), Center(), Center())
times = N.times

hmN = heatmap!(axN, times ./ year, zc, N[1, 1, 1:grid.Nz, 1:end]', interpolate = true, colormap = Reverse(:batlow))

hmP = heatmap!(axP, times ./ year, zc, P[1, 1, 1:grid.Nz, 1:end]', interpolate = true, colormap = Reverse(:batlow))

lines!(axSed, times ./ year, sed[1, 1, 1, :])

Colorbar(fig[1, 2], hmN, label = "Nutrient (mmol N / m³)")
Colorbar(fig[2, 2], hmP, label = "Phytoplankton (mmol N / m³)")

fig
```
</details>

We can see in this that some phytoplankton sink to the bottom, and are both remineralized back into nutrients and stored in the sediment.

### Final notes
When implementing a new model we recommend following a testing process as we have here, starting with a box model, then a column, and finally using it in a realistic physics scenarios. We have found this very helpful for spotting bugs that were proving difficult to decipher in other situations. You can also add `Individuals`, light attenuation models, and sediment models in a similar fashion.