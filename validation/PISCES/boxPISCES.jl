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

grid = RectilinearGrid( topology = (Flat, Flat, Flat), size = (), z = -100)
clock = Clock(time = zero(grid))

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR_func(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
MLD(t) = - 100
EU(t) = - 50
 # specify the nominal depth of the box for the PAR profile
# Modify the PAR based on the nominal depth and exponential decay
#PAR_func(t) = 18.0 # Modify the PAR based on the nominal depth and exponential decay
#PAR_func(t) = 18.0
PAR_func1(t) = PAR_func(t)/3
PAR_func2(t) = PAR_func(t)/3
PAR_func3(t)= PAR_func(t)/3

PAR = FunctionField{Center, Center, Center}(PAR_func, grid; clock)
PAR¹ = FunctionField{Center, Center, Center}(PAR_func1, grid; clock)
PAR² = FunctionField{Center, Center, Center}(PAR_func2, grid; clock)
PAR³ = FunctionField{Center, Center, Center}(PAR_func3, grid; clock)
zₘₓₗ = FunctionField{Center, Center, Center}(MLD, grid; clock)
zₑᵤ = FunctionField{Center, Center, Center}(EU, grid; clock)
nothing #hide

# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
model = BoxModel(; biogeochemistry = PISCES(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation((; PAR, PAR¹, PAR², PAR³)), mixed_layer_depth = zₘₓₗ, euphotic_layer_depth = zₑᵤ),
                   clock)

set!(model, P = 6.95, D = 6.95, Z = 0.695,  M = 0.695, Pᶜʰˡ = 1.671,  Dᶜʰˡ = 1.671, Pᶠᵉ =7e-6 * 1e9 / 1e6 * 6.95, Dᶠᵉ = 7e-6 * 1e9 / 1e6 * 6.95, Dˢⁱ = 1.162767, DOC = 0.0, POC = 0.0, GOC = 0.0, SFe = 7e-6 * 1e9 / 1e6 *1.256, BFe =7e-6 * 1e9 / 1e6 *1.256, NO₃ = 6.202, NH₄ = 0.25*6.202, PO₄ = 0.8722, Fe = 1.256, Si = 7.313, CaCO₃ = 0.29679635764590534, DIC = 2139.0, Alk = 2366.0, O₂ = 237.0, T = 14.0) #Using Copernicus Data (26.665, 14.)

simulation = Simulation(model; Δt = 90minutes, stop_time =5years)

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box.jld2", schedule = TimeInterval(10day), overwrite_existing = true)

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

fig = Figure(size = (2400, 3600), fontsize = 24)

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/4), Int(idx%4)], ylabel = "$name", xlabel = "years", xticks=(0:40)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end
push!(axs, Axis(fig[6,1], ylabel = "PAR", xlabel = "years", xticks=(0:40)))
lines!(axs[end], (0:10day:5years) / year, x -> PAR_func(x * year), linewidth = 3)
fi = length(timeseries.P)

#println("P = $(timeseries.P[fi]), D = $(timeseries.D[fi]), Z = $(timeseries.Z[fi]),  M = $(timeseries.M[fi]), Pᶜʰˡ = $(timeseries.Pᶜʰˡ[fi]),  Dᶜʰˡ = $(timeseries.Dᶜʰˡ[fi]), Pᶠᵉ = $(timeseries.Pᶠᵉ[fi]), Dᶠᵉ = $(timeseries.Dᶠᵉ[fi]), Dˢⁱ = $(timeseries.Dˢⁱ[fi]), DOC = $(timeseries.DOC[fi]), POC = $(timeseries.POC[fi]), GOC = $(timeseries.GOC[fi]), SFe = $(timeseries.SFe[fi]), BFe = $(timeseries.BFe[fi]), PSi = $(timeseries.PSi[fi]), NO₃ = $(timeseries.NO₃[fi]), NH₄ = $(timeseries.NH₄[fi]), PO₄ = $(timeseries.PO₄[fi]), Fe = $(timeseries.Fe[fi]), Si = $(timeseries.Si[fi]), CaCO₃ = $(timeseries.CaCO₃[fi]), DIC = $(timeseries.DIC[fi]), Alk = $(timeseries.Alk[fi]), O₂ = $(timeseries.O₂[fi]), T = 14.0")

Carbon_at_start = timeseries.P[1] + timeseries.D[1] + timeseries.Z[1] + timeseries.M[1] + timeseries.DOC[1] + timeseries.POC[1] + timeseries.GOC[1] + timeseries.DIC[1] + timeseries.CaCO₃[1]
Carbon_at_end = timeseries.P[fi] + timeseries.D[fi] + timeseries.Z[fi] + timeseries.M[fi] + timeseries.DOC[fi] + timeseries.POC[fi] + timeseries.GOC[fi] + timeseries.DIC[fi] + timeseries.CaCO₃[fi]
Iron_at_start = 10e-3*timeseries.Z[1] + 10e-3*timeseries.M[1] + timeseries.Pᶠᵉ[1] + timeseries.Dᶠᵉ[1] + timeseries.Fe[1] + timeseries.BFe[1] + timeseries.SFe[1]
Iron_at_end = 10e-3*timeseries.Z[fi] + 10e-3*timeseries.M[fi] + timeseries.Pᶠᵉ[fi] + timeseries.Dᶠᵉ[fi] + timeseries.Fe[fi] + timeseries.BFe[fi] + timeseries.SFe[fi]
Silicon_at_start = timeseries.Dˢⁱ[1] + timeseries.Si[1] + timeseries.PSi[1]
Silicon_at_End = timeseries.Dˢⁱ[fi] + timeseries.Si[fi] + timeseries.PSi[fi]
Nitrogen_at_start = 16/122*(timeseries.P[1] + timeseries.D[1] + timeseries.Z[1] + timeseries.M[1] + timeseries.DOC[1] + timeseries.POC[1] + timeseries.GOC[1]) + timeseries.NO₃[1] + timeseries.NH₄[1]
Nitrogen_at_end =  16/122*(timeseries.P[fi] + timeseries.D[fi] + timeseries.Z[fi] + timeseries.M[fi] + timeseries.DOC[fi] + timeseries.POC[fi] + timeseries.GOC[fi]) + timeseries.NO₃[fi] + timeseries.NH₄[fi]
Phosphates_at_start = 1/122*(timeseries.P[1] + timeseries.D[1] + timeseries.Z[1] + timeseries.M[1] + timeseries.DOC[1] + timeseries.POC[1] + timeseries.GOC[1]) + timeseries.PO₄[1]
Phosphates_at_end = 1/122*(timeseries.P[fi] + timeseries.D[fi] + timeseries.Z[fi] + timeseries.M[fi] + timeseries.DOC[fi] + timeseries.POC[fi] + timeseries.GOC[fi]) + timeseries.PO₄[fi]

println("Carbon at start = ", Carbon_at_start, " Carbon at end = ", Carbon_at_end)
println("Iron at start = ", Iron_at_start, " Iron at end = ", Iron_at_end)
println("Silicon at start = ", Silicon_at_start, " Silicon at end = ", Silicon_at_End)
println("Phosphates at start = ", Phosphates_at_start, " Phosphates at end = ", Phosphates_at_end)
println("Nitrogen at start = ", Nitrogen_at_start, " Nitrogen at end = ", Nitrogen_at_end)
fig