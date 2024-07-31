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
using Oceananigans: UpdateStateCallsite, TendencyCallsite

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
model = BoxModel(; biogeochemistry = PISCES(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR¹, PAR², PAR³)), flux_feeding_rate = 2.0e-3),
                   clock)

set!(model, NO₃ = 4.0, NH₄ = 0.1, P = 4.26, D = 4.26, Z = 0.426, M = 0.426, Pᶠᵉ = 3.50, Dᶠᵉ = 3.50, Pᶜʰˡ = 1.60, Dᶜʰˡ = 1.60, Dˢⁱ = 0.525, Fe = 0.00082410, O₂ = 264.0, Si = 4.557, Alk = 2360.0, PO₄ = 0.8114, DIC = 2000.0, CaCO₃ = 1.0, T = 14.0)
#set!(model, P = 2.7591635623893302, D = 0.11455452311786936, Z = 1.1582114117772278,  M = 4.063075192668189, Pᶜʰˡ = 0.2982091168733439,  Dᶜʰˡ = 0.016732978567852125, Pᶠᵉ = 16.728479580167, Dᶠᵉ = 4.546624294394084, Dˢⁱ = 0.02069887213660767, DOC = 0.006064167787693882, POC = 0.03005527112604406, GOC = 0.002767910547973816, SFe = 10.923656390454695, BFe = 0.3667994320927491, PSi = 0.8763378049500037, NO₃ = 2.6202550254764647, NH₄ = 1.8990251077677787, PO₄ = 0.8320698170409203, Fe = 6.668323320599245, Si = 4.493774109241143, CaCO₃ = 1.6345832401588156, DIC = 2001.88713443883, Alk = 2360.593403431815, O₂ = 263.90313885636107, T = 14.0)

simulation = Simulation(model; Δt = 5minutes, stop_time = 10days)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

function non_zero_fields!(model) 
    #for tracer in model.fields
      #  if tracer != model.fields[:T]
       #     tracer = tracer[1,1,1]
        #    getproperty(model, tracer) = max(0.0, getproperty(model, tracer))
        #end
    @inbounds for field in model.fields
        field[1, 1, 1] = max(0, field[1, 1, 1])
    end
    #end
    return nothing

end

simulation.callbacks[:progress] = Callback(prog, IterationInterval(10))
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