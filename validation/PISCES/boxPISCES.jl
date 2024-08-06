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
using Oceananigans.Fields: FunctionField

const year = years = 365day
nothing #hide

grid = BoxModelGrid()
clock = Clock(time = zero(grid))

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

const z = - 10 # specify the nominal depth of the box for the PAR profile
# Modify the PAR based on the nominal depth and exponential decay
PAR_func(t) = 300.0 # Modify the PAR based on the nominal depth and exponential decay

PAR_func1(t) = 100.0
PAR_func2(t) = 100.0
PAR_func3(t)= 100.0

PAR = FunctionField{Center, Center, Center}(PAR_func, grid; clock)
PAR¹ = FunctionField{Center, Center, Center}(PAR_func1, grid; clock)
PAR² = FunctionField{Center, Center, Center}(PAR_func2, grid; clock)
PAR³ = FunctionField{Center, Center, Center}(PAR_func3, grid; clock)

nothing #hide

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(; biogeochemistry = PISCES(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR¹, PAR², PAR³))),
                   clock)

#set!(model, NO₃ = 4.0, NH₄ = 0.1, P = 4.26, D = .426, Z = 0.426, M = 0.426, Pᶠᵉ = 3.50, Dᶠᵉ = 3.50, Pᶜʰˡ = .1, Dᶜʰˡ = .01, Dˢⁱ = 0.525, Fe = 0.00082410, O₂ = 264.0, Si = 4.557, Alk = 2360.0, PO₄ = 0.8114, DIC = 2000.0, CaCO₃ = 0.0001, T = 14.0)
set!(model, P = 3.4969507431335574, D = 3.045411522412732, Z = 0.06632989062825102,  M = 0.2013934196815033, Pᶜʰˡ = 0.1649236496634812,  Dᶜʰˡ = 0.1744075407918407, Pᶠᵉ = 0.629885753182277, Dᶠᵉ = 0.5811469701875243, Dˢⁱ = 0.6798967635852701, DOC = 0.00038765439868301653, POC = 0.0015526645806093746, GOC = 15.59902032438512, SFe = 0.0003450130201159402, BFe = 3.7215124823095893, PSi = 0.5404372397384611, NO₃ = 2.828750058046986, NH₄ = 0.5246932481313445, PO₄ = 0.7647402066361438, Fe = 2.240993210202666e-5, Si = 3.7394628646278316, CaCO₃ = 0.00933485588389991, DIC = 1994.2982703537236, Alk = 2361.5639068719747, O₂ = 272.5482524012977, T = 14.0)

simulation = Simulation(model; Δt = 5minutes, stop_time = 7hours)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(5minutes), overwrite_existing = true)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

function non_zero_fields!(model) 
    #for tracer in model.fields
      #  if tracer != model.fields[:T]
       #     tracer = tracer[1,1,1]
        #    getproperty(model, tracer) = max(0.0, getproperty(model, tracer))
        #end
    @inbounds for (idx, field) in enumerate(model.fields)
        if isnan(field[1,1,1])
            throw("$(keys(model.fields)[idx]) has gone NaN")
        else
            field[1, 1, 1] = max(0, field[1, 1, 1])
        end
        
    end
    #end
    return nothing

end

simulation.callbacks[:progress] = Callback(prog, IterationInterval(100))
simulation.callbacks[:non_zero_fields] = Callback(non_zero_fields!, callsite = UpdateStateCallsite())


@info "Running the model..."
run!(simulation)

# ## Load the output

times = FieldTimeSeries("box.jld2", "P").times

timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("box.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 7200), fontsize = 24)

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end

fi = length(timeseries.P)

println("P = $(timeseries.P[fi]), D = $(timeseries.D[fi]), Z = $(timeseries.Z[fi]),  M = $(timeseries.M[fi]), Pᶜʰˡ = $(timeseries.Pᶜʰˡ[fi]),  Dᶜʰˡ = $(timeseries.Dᶜʰˡ[fi]), Pᶠᵉ = $(timeseries.Pᶠᵉ[fi]), Dᶠᵉ = $(timeseries.Dᶠᵉ[fi]), Dˢⁱ = $(timeseries.Dˢⁱ[fi]), DOC = $(timeseries.DOC[fi]), POC = $(timeseries.POC[fi]), GOC = $(timeseries.GOC[fi]), SFe = $(timeseries.SFe[fi]), BFe = $(timeseries.BFe[fi]), PSi = $(timeseries.PSi[fi]), NO₃ = $(timeseries.NO₃[fi]), NH₄ = $(timeseries.NH₄[fi]), PO₄ = $(timeseries.PO₄[fi]), Fe = $(timeseries.Fe[fi]), Si = $(timeseries.Si[fi]), CaCO₃ = $(timeseries.CaCO₃[fi]), DIC = $(timeseries.DIC[fi]), Alk = $(timeseries.Alk[fi]), O₂ = $(timeseries.O₂[fi]), T = 14.0")

fig