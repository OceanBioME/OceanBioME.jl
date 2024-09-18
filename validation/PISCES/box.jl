# # [Box model](@id box_example)
# In this example we setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME, Oceananigans, Oceananigans.Units

using Oceananigans.Fields: FunctionField, ConstantField

using OceanBioME.Models.PISCESModel: SimpleIron, NitrateAmmonia

const year = years = 365day

grid = RectilinearGrid( topology = (Flat, Flat, Flat), size = (), z = -10)
clock = Clock(time = zero(grid))

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
@inline surface_PAR(t) = 300 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 10
@inline temp(t) = 2.4 * (1 + cos(t * 2π / year + 50days)) + 8

@inline PAR_component(t) = surface_PAR(t) * exp(-10/2) / 3

PAR₁ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)
PAR₂ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)
PAR₃ = FunctionField{Nothing, Nothing, Center}(PAR_component, grid; clock)

PAR = PAR₁ + PAR₂ + PAR₃

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR₁, PAR₂, PAR₃))

mixed_layer_depth = ConstantField(-100)
euphotic_depth = ConstantField(-100)
mean_mixed_layer_vertical_diffusivity = ConstantField(1e-2)
mean_mixed_layer_light = PAR
yearly_maximum_silicate = ConstantField(0) # turn off the contribution from enhanced requirments

biogeochemistry = PISCES(; grid, 
                           sinking_speeds = (POC = 0, GOC = 0), 
                           light_attenuation, 
                           mixed_layer_depth, 
                           euphotic_depth,
                           yearly_maximum_silicate,
                           mean_mixed_layer_light,
                           mean_mixed_layer_vertical_diffusivity,
                           iron = SimpleIron(0),
                           nitrogen = NitrateAmmonia(maximum_fixation_rate = 0.0))

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(; grid, biogeochemistry, clock, prescribed_tracers = (; T = temp))

set!(model, P = 6.95, D = 6.95, Z = 0.695,  M = 0.695, 
            PChl = 1.671,  DChl = 1.671, 
            PFe = 7e-6 * 1e9 / 1e6 * 6.95, DFe = 7e-6 * 1e9 / 1e6 * 6.95, 
            DSi = 1.162767, 
            NO₃ = 6.202, NH₄ = 0.25*6.202, 
            PO₄ = 0.8722, Fe = 1.256, Si = 7.313, 
            CaCO₃ = 100,
            DIC = 2139.0, Alk = 2366.0, 
            O₂ = 237.0, S = 35, T = 10) 

            
simulation = Simulation(model; Δt = 20minutes, stop_time = 4years)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(0.5day), overwrite_existing = true)

PAR_field = Field(biogeochemistry.light_attenuation.fields[1])
simulation.output_writers[:par] = JLD2OutputWriter(model, (; PAR = PAR_field); filename = "box_light.jld2", schedule = TimeInterval(0.5day), overwrite_existing = true)


prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

#NaN checker function, could be removed if confident of model stability
function non_zero_fields!(model) 
    @inbounds for (idx, field) in enumerate(model.fields)
        if isnan(field[1,1,1])
            throw("$(keys(model.fields)[idx]) has gone NaN")
        else
            field[1, 1, 1] = max(0, field[1, 1, 1])
        end
        
    end
    return nothing
end

simulation.callbacks[:progress] = Callback(prog, TimeInterval(182days))
#simulation.callbacks[:non_zero_fields] = Callback(non_zero_fields!, callsite = UpdateStateCallsite())

@info "Running the model..."
run!(simulation)

# ## Load the output
timeseries = FieldDataset("box.jld2")

times = timeseries.fields["P"].times

PAR_timeseries = FieldTimeSeries("box_light.jld2", "PAR")
# ## And plot
using CairoMakie

fig = Figure(size = (2400, 3600), fontsize = 24)

axs = []
for name in Oceananigans.Biogeochemistry.required_biogeochemical_tracers(biogeochemistry)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/4), Int(idx%4)], ylabel = "$name", xlabel = "years", xticks=(0:40)))
    lines!(axs[end], times[731:end] / year, timeseries["$name"][731:end], linewidth = 3)
end

#Plotting the function of PAR
push!(axs, Axis(fig[6, 2], ylabel = "PAR", xlabel = "years", xticks=(0:40)))
lines!(axs[end], times[731:end]/year, PAR_timeseries[731:end], linewidth = 3)

fig