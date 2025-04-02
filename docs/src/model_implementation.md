# [Implementing new models](@id model_implementation)

Here we describe how OceanBioME defines biogeochemical (BGC) models, how this varies from Oceananigans, and how to implement your own model.

## Model structure
OceanBioME BGC models are `struct`s of type `ContinuousFormBiogeochemistry`, which is of abstract type `AbstractContinuousFormBiogeochemistry` from Oceananigans. In Oceananigans this describes BGC models which are defined using continuous functions (depending continuously on ``x``, ``y``, and ``z``) rather than discrete functions (depending on ``i``, ``j``, ``k``). This allows the user to implement the BGC model equations without worrying about details of the grid or discretization, and then Oceananigans handles the rest.

OceanBioME's `ContinuousFormBiogeochemistry` adds a layer on top of this which makes it easy to add [light attenuation models](@ref light), [sediment](@ref sediment), and [biologically active particles](@ref individuals) (or individual-based models). OceanBioME's `ContinuousFormBiogeochemistry` includes parameters in which the types of these components are stored. This means that these model components will automatically be integrated into the BGC model without having to add new methods to call Oceananigans functions. 

## Implementing a model

The nature of multiple dispatch in Julia means that we define new BGC models as new types. You can then define [methods](https://docs.julialang.org/en/v1/manual/methods/) to this type which are used by OceanBioME and Oceananigans to integrate the model.

### The basics

For this example we are going to implement the simple Nutrient-Phytoplankton model similar to that used in [Chen2015](@citep), although we neglect the nutrient in/outflow terms since they may be added as [boundary conditions](https://clima.github.io/OceananigansDocumentation/stable/model_setup/boundary_conditions/), and modified to conserve nitrogen.

The first step is to import the abstract type from OceanBioME, some units from Oceananigans (for ease of parameter definition), and [`import`](https://docs.julialang.org/en/v1/manual/faq/#What-is-the-difference-between-%22using%22-and-%22import%22?) some functions from Oceananigans in order to add methods to:

```@example implementing
using OceanBioME, Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity
```

We then define our `struct` with the model parameters, as well as slots for the particles, light attenuation, and sediment models:

```@example implementing
@kwdef struct NutrientPhytoplankton{FT, W} <: AbstractContinuousFormBiogeochemistry
            base_growth_rate :: FT = 1.27 / day              # 1 / seconds
    nutrient_half_saturation :: FT = 0.025 * 1000 / 14       # mmol N / m³
       light_half_saturation :: FT = 300.0                   # micro einstein / m² / s
        temperature_exponent :: FT = 0.24                    # 1
     temperature_coefficient :: FT = 1.57                    # 1
         optimal_temperature :: FT = 28.0                    # °C
              mortality_rate :: FT = 0.15 / day              # 1 / seconds
     crowding_mortality_rate :: FT = 0.004 / day / 1000 * 14 # 1 / seconds / mmol N / m³
            sinking_velocity :: W  = 2 / day
end
```

Here, we use descriptive names for the parameters. Below, each of these parameters correspond to a symbol (or letter) which is more convenient mathematically and when defining the BGC model functions. In the above code we used `@kwdef` to set default values for the models so that we don't have to set all of these parameters each time we use the model. The default parameter values can optionally be over-ridden by the user when running the model. We have also included a `sinking_velocity` field in the parameter set to demonstrate how we can get tracers (e.g. detritus) to sink. We also need to define some functions so that OceanBioME and Oceananigans know what tracers and auxiliary fields (e.g. light intensity) we use:

```@example implementing
required_biogeochemical_tracers(::NutrientPhytoplankton) = (:N, :P, :T)

required_biogeochemical_auxiliary_fields(::NutrientPhytoplankton) = (:PAR, )
```

Next, we define the functions that specify how the phytoplankton ``P`` evolve. In the absence of advection and diffusion (both of which are handled by Oceananigans), we want the phytoplankton to evolve at the rate given by:
```math
\frac{\partial P}{\partial t} = \mu g(T) f(N) h(PAR) P - mP - bP^2,
```

where ``\mu`` corresponds to the parameter `base_growth_rate`, ``m`` corresponds to the parameter `mortality_rate`, and ``b`` corresponds to the parameter `crowding_mortality_rate`. Here, the functions ``g``, ``f``, and ``h`` are defined by:
```math
\begin{align}
g(T) &= c_1\exp\left(-c_2|T - T_{opt}|\right),\\
f(N) &= \frac{N}{k_N + N},\\
h(PAR) &= \frac{PAR}{k_P + PAR},
\end{align}
```
where ``c_1`` corresponds to `temperature_coefficient`,  ``c_2`` corresponds to `temperature_exponent`, ``T_{opt}`` corresponds to `optimal_temperature`, ``k_N`` corresponds to `nutrient_half_saturation`, and ``k_P`` corresponds to `light_half_saturation`. 

We turn this into a function for our model by writing:

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

The first parameter `::Val{:P}` is a special [value type](http://www.jlhub.com/julia/manual/en/function/Val) that allows this function to be dispatched when it is given the value `Val(:P)`. This is how Oceananigans tells the model which forcing function to use. At the start of the `NutrientPhytoplankton` function we unpack some parameters from the model, then calculate each term, and return the total change (the gain minus the loss). 

For this model, the nutrient evolution can be inferred from the rate of change of phytoplankton. Since this is a simple two variable model and the total concentration is conserved, 
```math
\frac{\partial N}{\partial t} = - \frac{\partial P}{\partial t}
```
Hence, we define the nutrient forcing using as the negative of the phytoplankton forcing
```@example implementing
@inline (bgc::NutrientPhytoplankton)(::Val{:N}, args...) = -bgc(Val(:P), args...)
```

Now we can run an example similar to the [LOBSTER box model example](@ref box_example):

```@example implementing
using OceanBioME, Oceananigans.Units
using Oceananigans.Fields: FunctionField

const year = years = 365days

@inline PAR⁰(t) = 500 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

clock = Clock(; time = 0.0)

z = -10 # specify the nominal depth of the box for the PAR profile
@inline PAR_func(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay 

PAR = FunctionField{Center, Center, Center}(PAR_func, BoxModelGrid(); clock)

@inline temp(t) = 2.4 * cos(t * 2π / year + 50days) + 26

biogeochemistry = Biogeochemistry(NutrientPhytoplankton(); 
                                  light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PAR))

model = BoxModel(; biogeochemistry,
                   prescribed_tracers = (; T = temp),
                   clock)

set!(model, N = 15, P = 15)

simulation = Simulation(model; Δt = 5minutes, stop_time = 5years)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box_np.jld2", schedule = TimeInterval(10days), overwrite_existing = true)

# ## Run the model (should only take a few seconds)
@info "Running the model..."
run!(simulation)
```

We can then visualise this:

```@example implementing
P = FieldTimeSeries("box_np.jld2", "P")
N = FieldTimeSeries("box_np.jld2", "N")

times = P.times

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 480), fontsize = 20)

axN= Axis(fig[1, 1], ylabel = "Nutrient \n(mmol N / m³)")
lines!(axN, times / year, N[1, 1, 1, :], linewidth = 3)

axP = Axis(fig[1, 2], ylabel = "Phytoplankton \n(mmol N / m³)")
lines!(axP, times / year, P[1, 1, 1, :], linewidth = 3)

axPAR= Axis(fig[2, 1], ylabel = "PAR (einstein / m² / s)", xlabel = "Time (years)")
lines!(axPAR, times / year, PAR_func.(times), linewidth = 3)

axT = Axis(fig[2, 2], ylabel = "Temperature (°C)", xlabel = "Time (years)")
lines!(axT, times / year, temp.(times), linewidth = 3)

fig
```

So now we know it works.

### Phytoplankton sinking

Now that we have a fully working BGC model we might want to add some more features. Another aspect that is easy to add is negative buoyancy (sinking). To-do this all we do is add a method to the Oceananigans function `biogeochemical_drift_velocity`, and we use `::Val{:P}` to specify that only phytoplankton will sink. Above, we set the default value of the parameter `bgc.sinking_velocity`. We can override this when we call the BGC model like `NutrientPhytoplankton(; light_attenuation_model, sinking_velocity = 1/day)`. Note that before using `biogeochemical_drift_velocity`, we need to import several `Fields` from Oceananigans:

```@example implementing
using Oceananigans.Fields: ZeroField, ConstantField

biogeochemical_drift_velocity(bgc::NutrientPhytoplankton, ::Val{:P}) = 
    (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocity)
```

### Sediment model coupling

Another aspect that OceanBioME includes is sediment models. Doing this varies between sediment models, but for the most generic and simplest, all we need to do is add methods to two functions:

```@example implementing
using OceanBioME.Sediments: sinking_flux

import OceanBioME.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

@inline nitrogen_flux(i, j, k, grid, advection, bgc::NutrientPhytoplankton, tracers) =
     sinking_flux(i, j, k, grid, advection, Val(:P), bgc, tracers)
                 
@inline carbon_flux(i, j, k, grid, advection, bgc::NutrientPhytoplankton, tracers) = nitrogen_flux(i, j, k, grid, advection, bgc, tracers) * 6.56

@inline remineralisation_receiver(::NutrientPhytoplankton) = :N

@inline sinking_tracers(::NutrientPhytoplankton) = (:P, )
```

### Putting it together

Now that we have added these elements we can put it together into another simple example:
```@example implementing
using Oceananigans, OceanBioME
using OceanBioME.Sediments: InstantRemineralisation

# define some simple forcing

@inline surface_PAR(t) = 200 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline ∂ₜT(z, t) = - 2π / year * sin(t * 2π / year + 50days)

@inline κₚ(z) = 1e-2 * (1 + tanh((z - 50) / 10)) / 2 + 1e-4

# define the grid

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (32, ), x = 1, y = 1, z = (-100, 0))

# setup the biogeochemical model

light_attenuation = TwoBandPhotosyntheticallyActiveRadiation(; grid, surface_PAR)

sediment = InstantRemineralisation(; grid)

sinking_velocity = ZFaceField(grid)

w_sink(z) = 2 / day * tanh(z / 5)

set!(sinking_velocity, w_sink)

negative_tracer_scaling = ScaleNegativeTracers((:N, :P))

biogeochemistry = Biogeochemistry(NutrientPhytoplankton(; sinking_velocity);
                                  light_attenuation,
                                  sediment,
                                  modifiers = negative_tracer_scaling) 

κ = CenterField(grid)

set!(κ, κₚ)

# put the model together

model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              closure = ScalarDiffusivity(ν = κ; κ), 
                              forcing = (; T = ∂ₜT))

set!(model, P = 0.01, N = 15, T = 28)

# run

simulation = Simulation(model, Δt = 9minutes, stop_time = 1years)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "column_np.jld2",
                                                       schedule = TimeInterval(1day),
                                                       overwrite_existing = true)

simulation.output_writers[:sediment] = JLD2OutputWriter(model, model.biogeochemistry.sediment.fields,
                                                        indices = (:, :, 1),
                                                        filename = "column_np_sediment.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)

run!(simulation)
```

We can then visualise this:

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

hmN = heatmap!(axN, times ./ year, zc, N[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

hmP = heatmap!(axP, times ./ year, zc, P[1, 1, 1:grid.Nz, 1:end]',
               interpolate = true, colormap = Reverse(:batlow))

lines!(axSed, times ./ year, sed[1, 1, 1, :])

Colorbar(fig[1, 2], hmN, label = "Nutrient (mmol N / m³)")
Colorbar(fig[2, 2], hmP, label = "Phytoplankton (mmol N / m³)")

fig
```

We can see in this that some phytoplankton sink to the bottom, and are both remineralized back into nutrients and stored in the sediment.

### Running on a GPU

In order to run a BGC model on a GPU, the BGC model must first be `adapted`. After the definition of the BGC `struct`, we need to write:

```
using Adapt

import Adapt: adapt_structure

Adapt.adapt_structure(to, bgc::NutrientPhytoplankton) = NutrientPhytoplankton(adapt(to, bgc.base_growth_rate),
                                                                              adapt(to, bgc.nutrient_half_saturation),
                                                                              adapt(to, bgc.light_half_saturation),
                                                                              adapt(to, bgc.temperature_exponent),
                                                                              adapt(to, bgc.temperature_coefficient),
                                                                              adapt(to, bgc.optimal_temperature),
                                                                              adapt(to, bgc.mortality_rate),
                                                                              adapt(to, bgc.crowding_mortality_rate),
                                                                              adapt(to, bgc.sinking_velocity))
```

Also, in order for `ScaleNegativeTracers` to work on a GPU, we must add `grid` as one of the input arguments. We replace the definition of `negative_tracer_scaling` with:
```
negative_tracer_scaling = ScaleNegativeTracers((:N, :P), grid)
```
We can then add `GPU()` to the definition of `grid` in the usual way, and the column model above will be able to run on a GPU. 

### Final notes
When implementing a new model we recommend following a testing process as we have here, starting with a box model, then a column, and finally using it in a realistic physics scenarios. We have found this very helpful for spotting bugs that were proving difficult to decipher in other situations. You can also add `Individuals`, light attenuation models, and sediment models in a similar fashion.
